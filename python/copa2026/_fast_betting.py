"""Batched history-only candidate calculations for the normalized strategy.

The mathematical formulas and cautious rule are unchanged. For streams longer
than 5000 positions the reference rolling-window implementation is retained.
"""
from __future__ import annotations
import numpy as np
from .betting import (_probability_array,prepare_knn,evaluate_knn,_combination,
    CANDIDATE_NAMES,STATUS_NAMES,_HIST_CENTERS)


def _tail_z(sorted_history,k,counts=None):
    b,n=sorted_history.shape
    counts=np.full(b,n,dtype=np.int64) if counts is None else np.asarray(counts,dtype=np.int64)
    coefficient=(k-1)/(2*counts)
    first=sorted_history[:,:n-k+1]
    last=sorted_history[:,k-1:]
    valid=np.arange(n-k+1)[None,:] < (counts-k+1)[:,None]
    with np.errstate(over='ignore',divide='ignore',invalid='ignore',under='ignore'):
        transitions=sorted_history[:,:n-k]+(sorted_history[:,k:]-sorted_history[:,:n-k])/2
        divergent=np.any(valid&(first==last),axis=1)
        minimum_radius=np.min(np.where(valid,last-first,np.inf),axis=1)/2
        maximum_raw=coefficient/minimum_radius
        early_failure=(minimum_radius<=0)|~np.isfinite(maximum_raw)
        left=np.concatenate((np.zeros((b,1)),transitions),axis=1)
        right=np.concatenate((transitions,np.ones((b,1))),axis=1)
        right[np.arange(b),counts-k]=1.
        middle=first+(last-first)/2
        low_end=np.minimum(right,middle)
        low_width=np.where(valid,np.maximum(0,low_end-left),0.)
        low_radius=last-low_end
        high_start=np.maximum(left,middle)
        high_width=np.where(valid,np.maximum(0,right-high_start),0.)
        high_radius=high_start-first
        low_active=low_width>0;high_active=high_width>0
        low_ratio=np.zeros_like(low_width);high_ratio=np.zeros_like(high_width)
        np.divide(low_width,low_radius,out=low_ratio,where=low_active)
        np.divide(high_width,high_radius,out=high_ratio,where=high_active)
        low_terms=np.log1p(low_ratio);high_terms=np.log1p(high_ratio)
        z=np.array([coef*np.sum(np.concatenate((lo[la],hi[ha])))
                    for coef,lo,la,hi,ha in zip(coefficient,low_terms,low_active,high_terms,high_active)])
        nonfinite_terms=np.any(~np.isfinite(low_terms),axis=1)|np.any(~np.isfinite(high_terms),axis=1)
        late_failure=~np.isfinite(z)|(z<=0)|nonfinite_terms|~np.isfinite(maximum_raw/z)
    numerical=~divergent&(early_failure|late_failure)
    fallback=divergent|numerical
    z[divergent]=np.inf
    z[~divergent&early_failure]=np.nan
    status=np.ones(b,dtype=np.uint8);status[divergent]=3;status[numerical]=4
    return z,fallback,status


def _histogram_candidates(pvalues):
    n=len(pvalues);row=np.arange(n);start=np.maximum(0,row-200)
    all_counts=[];occupied=[]
    for bins,centers in enumerate(_HIST_CENTERS,1):
        # Match HIST's left + diff/2 edges, eps shift and HISTC inclusion.
        edges=centers[:-1]+np.diff(centers)/2
        edges=edges+np.spacing(edges)
        membership=np.searchsorted(edges,pvalues,side='right')
        one_hot=np.zeros((n,bins),dtype=np.int64)
        one_hot[row,membership]=1
        prefix=np.concatenate((np.zeros((1,bins),dtype=np.int64),np.cumsum(one_hot,axis=0)),axis=0)
        counts=prefix[row]-prefix[start]
        all_counts.append(counts)
        occupied.append(np.all(counts>0,axis=1))
    bets=np.ones((n,3));query=np.maximum(pvalues,0.00001)
    for col,maximum in enumerate((5,10,15)):
        selected_bins=np.ones(n,dtype=int)
        for bins in range(1,maximum+1):selected_bins[occupied[bins-1]]=bins
        for bins in np.unique(selected_bins[2:]):
            indices=np.flatnonzero((selected_bins==bins)&(row>=2))
            counts=all_counts[bins-1][indices]
            density=counts/np.sum(counts,axis=1,keepdims=True)*bins
            centers=_HIST_CENTERS[bins-1]
            distances=np.abs(centers[None,:]-query[indices,None])
            active=distances<1/bins
            weights=np.where(active,1/bins-distances,0)
            weights/=np.sum(weights,axis=1,keepdims=True)
            bets[indices,col]=np.sum(weights*density,axis=1)
    return bets


def fast_combination(pvalues,epsilon=100.,threshold=100.,return_diagnostics=False):
    pvalues=_probability_array(pvalues,vector=True).reshape(-1)
    for name,value in (('epsilon',epsilon),('threshold',threshold)):
        if not np.isscalar(value) or not np.isfinite(value) or value<=0:raise ValueError(f'{name} must be finite and positive')
    n=len(pvalues)
    if n<3 or n>5000:
        return _combination(pvalues,epsilon,threshold,return_diagnostics)
    values=np.ones((n,7));zs=np.full((n,3),np.nan);fallbacks=np.zeros((n,3),bool)
    statuses=np.zeros((n,3),np.uint8);raws=np.full((n,3),np.nan)
    for i in range(2,min(n,15)):
        history=pvalues[max(0,i-100):i]
        for j,k in enumerate((5,10,15)):
            model=prepare_knn(history,k);value,raw=evaluate_knn(model,pvalues[i])
            zs[i,j]=model.Z;fallbacks[i,j]=model.fallback;statuses[i,j]=STATUS_NAMES.index(model.status)
            raws[i,j]=raw;values[i,3-j]=value
    if n>15:
        indices=np.arange(15,n)
        counts=np.minimum(indices,100)
        history_indices=indices[:,None]-np.arange(100,0,-1)[None,:]
        history=np.where(history_indices>=0,pvalues[np.maximum(history_indices,0)],np.inf)
        sorted_history=np.sort(history,axis=1)
        # Padded infinities are masked; each real entry precedes its query.
        for j,k in enumerate((5,10,15)):
            zs[15:,j],fallbacks[15:,j],statuses[15:,j]=_tail_z(sorted_history,k,counts)
        # Match MATLAB Euclidean-distance squaring, including underflow.
        with np.errstate(under="ignore"):
            distances=np.sqrt(np.square(sorted_history-pvalues[15:,None]))
        distances=np.partition(distances,(4,9,14),axis=1)
        radius=distances[:,(4,9,14)]
        dist=2*radius;dist[dist==0]=1
        with np.errstate(over='ignore',divide='ignore',invalid='ignore'):
            raw=np.array([4.,9.,14.])[None,:]/(counts[:,None]*dist)
            normalized=np.divide(raw,zs[15:])
        normalized=np.where(fallbacks[15:],1.,normalized)
        if np.any(~np.isfinite(normalized)|(normalized<0)):raise FloatingPointError('Unrepresentable normalized factor')
        raws[15:]=raw;values[15:,1:4]=normalized[:,::-1]
    values[:,4:7]=_histogram_candidates(pvalues)
    with np.errstate(over='ignore',divide='ignore',invalid='ignore'):
        candidate_log=np.cumsum(np.log(values),axis=0)
        running_min=np.minimum.accumulate(candidate_log,axis=0)
        gains=np.zeros_like(candidate_log)
        gains[1:]=candidate_log[:-1]-running_min[:-1]
        best=np.argmax(gains,axis=1)
        selected=np.where(gains[np.arange(n),best]>np.log(epsilon),best+1,0)
        selected[:2]=0
        factors=np.where(selected>0,values[np.arange(n),best],1.)
        wealth=np.cumprod(factors)
        logs=np.cumsum(np.log(factors))
    diagnostics=None
    if return_diagnostics:
        diagnostics={'k':np.array([5,10,15]),'Z':zs,'fallback':fallbacks,'statusCode':statuses,
                     'statusNames':STATUS_NAMES,'selected':selected,'rawKnn':raws,
                     'candidateValues':values,'candidateLog':candidate_log,
                     'candidateNames':CANDIDATE_NAMES,'epsilon':epsilon,'W':5000,'threshold':threshold}
    return wealth,logs,diagnostics


  function [results, logResults, diagnostics] = combination(pvalues,epsilon,significance_l)
%COMBINATION Cautious betting with normalized kNN density candidates.
%   [RESULTS,LOGRESULTS,DIAGNOSTICS] = COMBINATION(PVALUES,EPSILON,THRESHOLD)
%   accepts a finite vector of p-values in [0,1] and returns row vectors of
%   linear capital and accumulated natural log capital, plus diagnostics.
%   EPSILON defaults to 100 and controls cautious activation; the third
%   argument defaults to capital threshold 100, not a significance level.
%   Alarm detection is performed by callers; this function does not stop.
%
%   The first two positions are startup/no-bet positions. kNN candidates use
%   up to 100 past p-values; histogram candidates use up to 200. Selection
%   uses past candidate log-capital drawups over a separate W=5000 window.
%   Histogram candidates use up to 5, 10 and 15 bins, reducing the bin count
%   when necessary to remove empty historical bins. Their formulas are kept.
%
%   Candidate columns are uniform, kNN15, kNN10, kNN5, hist5, hist10, hist15.
%   Selected code 0 means no bet. Normalizer/raw-kNN diagnostic columns use
%   k=[5 10 15], so their order differs from the candidate columns.
%   Normalizers and fallback decisions depend only on historical p-values.
%   All candidate histories are updated, including candidates not selected.
%
%   Requires KNN_BETTING_PREPARE and Statistics and Machine Learning Toolbox.
%   Analytical antiderivatives are evaluated in floating-point arithmetic.
%   Extremely tiny historical spans can make KNNSEARCH distances underflow
%   to zero, leaving its evaluation inconsistent with the analytic integral.
%   That known edge case is documented here; it is not corrected by comments.
%   LOGRESULTS avoids overflow/underflow of the linear product, but does not
%   repair errors in an individual betting factor. Normalization alone does
%   not establish all assumptions of an anytime-valid false-alarm guarantee.
if nargin < 2, epsilon = 100; end
if nargin < 3, significance_l = 100; end
validateattributes(pvalues, {'numeric'}, {'real','finite','>=',0,'<=',1});
if ~isvector(pvalues) && ~isempty(pvalues)
    error('combination:InvalidShape','pvalues must be a vector.');
end
validateattributes(epsilon, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(significance_l, {'numeric'}, {'scalar','real','finite','positive'});
pvalues = double(pvalues(:).');
n = numel(pvalues);
diagnostics = struct('k',[5 10 15], 'Z',nan(n,3), ...
    'fallback',false(n,3), 'statusCode',zeros(n,3,'uint8'), ...
    'statusNames',{{'startup','analytic','uniform_small_sample','uniform_divergent','uniform_numerical_failure'}}, ...
    'selected',zeros(n,1), 'rawKnn',nan(n,3), ...
    'candidateValues',ones(n,7), 'candidateLog',zeros(n,7), ...
    'candidateNames',{{'uniform','knn15','knn10','knn5','hist5','hist10','hist15'}}, ...
    'epsilon',epsilon, 'W',5000, 'threshold',significance_l);
% Accumulate log capital independently of the potentially overflowing linear product.
 logResults = zeros(1, numel(pvalues));
 if numel(pvalues) < 3
     results = ones(1, numel(pvalues));
     return;
 end
 choose=[];
  S5(1:length(pvalues))=1;S_5(1:length(pvalues))=1;
 
S2(1:length(pvalues))=1;S_2(1:length(pvalues))=1;
S3(1:length(pvalues))=1;
S4(1:length(pvalues))=1;S_4(1:length(pvalues))=1;


            S(1:length(pvalues))=1;
            g(1:length(pvalues))=1;
            f(1:length(pvalues))=1;
            f2(1:length(pvalues))=1;
            f3(1:length(pvalues))=1;
            h(1:length(pvalues))=1;
            h2(1:length(pvalues))=1;
            h3(1:length(pvalues))=1;
            h_(1:length(pvalues))=1;
            h_2(1:length(pvalues))=1;
            h_3(1:length(pvalues))=1;
            p=pvalues(1);
             
            S(1)=g(1);% Final selected martingale in linear capital.
            S2(1)=g(1);% kNN15 candidate; converted to log capital below.
            S3(1)=g(1);% Uniform candidate; converted to log capital below.
            S4(1)=g(1);
            S5(1)=g(1);
            S2=log(S2);
             S3=log(S3);
            S4=log(S4);
             S5=log(S5);
S_2(1)=g(1); S_2=log(S_2);S_5(1)=g(1); S_5=log(S_5);S_4(1)=g(1); S_4=log(S_4);


           for i=3:length(pvalues)
 pvalues_e=[  pvalues(max(1,i-100):i-1)  ]';% At most 100 historical p-values; current value excluded.

N=numel(pvalues_e);
% Prepare the three normalizers from history before evaluating the current p-value.
% Preserve the original nearest-neighbour definition and interval boundaries.
models = {knn_betting_prepare(pvalues_e,5), ...
          knn_betting_prepare(pvalues_e,10), ...
          knn_betting_prepare(pvalues_e,15)};
for j=1:3
    diagnostics.Z(i,j)=models{j}.Z;
    diagnostics.fallback(i,j)=models{j}.fallback;
    diagnostics.statusCode(i,j)=uint8(find(strcmp(models{j}.status,diagnostics.statusNames),1)-1);
end
% Uniform betting candidate: constant density g(i)=1.
 g(i)=1;  max(g(i),10^(-10));
 m=1;m2=1;m3=1;m4=1;m_2=1;m_3=1;m_4=1;
               counter=1;
            

               
               m=S3(i-1)-min(S3(max(i-5000,1):i-1));
               counter=counter+1;

 % Histogram candidate: up to 5 occupied bins, using at most 200 past values.
 a1=0;  
              no_bars_=5+1;
              while(~isempty(find(a1==0))&no_bars_>0)
               no_bars_=no_bars_-1;
              if (no_bars_>1)
               a1=hist(pvalues(max(1,i-200):i-1),[1/(no_bars_*2):1/no_bars_:1]);
              else
                a1=hist(pvalues(max(1,i-200):i-1),1);   
              end
              end
               centers=[1/(no_bars_*2):1/no_bars_:1]; 
              a1=a1/sum(a1)*length(a1);
                
                p=pvalues(i);
                p=max(p,0.00001);
                h_(i)=a1(ceil(p*length(a1)));
                
                if(~isempty(find(a1==0)))
                    h_(i)=1;
                    
                end               
               centers2=centers-p;
          tmp=find(abs(centers2)<1/no_bars_);
          centers=centers(tmp);
                          h_(i)=(1/no_bars_-abs(centers-p))/sum(1/no_bars_-(abs(centers-p)))*a1(tmp)';

                m_2=S_2(i-1)-min(S_2(max(i-5000,1):i-1));
% Histogram candidate: up to 10 occupied bins, with unchanged interpolation.
  a1=0;  
              no_bars_=10+1;
              while(~isempty(find(a1==0))&no_bars_>0)
               no_bars_=no_bars_-1;
              if (no_bars_>1)
               a1=hist(pvalues(max(1,i-200):i-1),[1/(no_bars_*2):1/no_bars_:1]);
              else
                a1=hist(pvalues(max(1,i-200):i-1),1);   
              end
              end
                 centers=[1/(no_bars_*2):1/no_bars_:1];
              a1=a1/sum(a1)*length(a1);
                
                p=pvalues(i);
                p=max(p,0.00001);
                
                if(~isempty(find(a1==0)))
                    h_2(i)=1;
                    
                end               
               centers2=centers-p;
          tmp=find(abs(centers2)<1/no_bars_);
          centers=centers(tmp);
                          h_2(i)=(1/no_bars_-abs(centers-p))/sum(1/no_bars_-(abs(centers-p)))*a1(tmp)';
                m_3=S_5(i-1)-min(S_5(max(i-5000,1):i-1));              
  
% Histogram candidate: up to 15 occupied bins, with unchanged interpolation.
  a1=0;  
              no_bars_=15+1;
              while(~isempty(find(a1==0))&no_bars_>0)
               no_bars_=no_bars_-1;
              if (no_bars_>1)
               a1=hist(pvalues(max(1,i-200):i-1),[1/(no_bars_*2):1/no_bars_:1]);
              else
                a1=hist(pvalues(max(1,i-200):i-1),1);   
              end
              end
                                          centers=[1/(no_bars_*2):1/no_bars_:1];
               a1=a1/sum(a1)*length(a1);
                
                p=pvalues(i);
                p=max(p,0.00001);
                
                if(~isempty(find(a1==0)))
                    h_3(i)=1;
                    
                end     
                centers2=centers-p;
          tmp=find(abs(centers2)<1/no_bars_);
          centers=centers(tmp);
                          h_3(i)=(1/no_bars_-abs(centers-p))/sum(1/no_bars_-(abs(centers-p)))*a1(tmp)';
                m_4=S_4(i-1)-min(S_4(max(i-5000,1):i-1));                
               % kNN15 candidate (candidate column 2; normalizer column 3).
 
               k=15;
   [mIdx,d] = knnsearch(pvalues_e,pvalues(i),'K',min(N,k))   ;        
dist=2*max(d);
if dist==0
dist=1; % Original zero-distance convention, also reached after distance underflow.
end
 
   h(i)=min(N-1,k-1)/(N*dist);
   diagnostics.rawKnn(i,3)=h(i);
   if models{3}.fallback
       h(i)=1; % Uniform fallback selected by the history-only preparation helper.
   else
       h(i)=h(i)/models{3}.Z;
   end
                
                              
               
                m2=S2(i-1)-min(S2(max(i-5000,1):i-1));
  % kNN10 candidate (candidate column 3; normalizer column 2).
k=10;
mIdx=mIdx(1:min(k,length(mIdx)));   
d =d(1:min(k,length(d)));        
dist=2*max(d);
if dist==0
dist=1; % Original zero-distance convention, also reached after distance underflow.
end
 
   h2(i)=min(N-1,k-1)/(N*dist);
   diagnostics.rawKnn(i,2)=h2(i);
   if models{2}.fallback
       h2(i)=1; % Uniform fallback selected by the history-only preparation helper.
   else
       h2(i)=h2(i)/models{2}.Z;
   end

               
               
                m3=S5(i-1)-min(S5(max(i-5000,1):i-1));              
  
% kNN5 candidate (candidate column 4; normalizer column 1).
k=5;
   mIdx=mIdx(1:min(k,length(mIdx)));   
d =d(1:min(k,length(d)));       
dist=2*max(d);
if dist==0
dist=1; % Original zero-distance convention, also reached after distance underflow.
end
 
   h3(i)=min(N-1,k-1)/(N*dist);
   diagnostics.rawKnn(i,1)=h3(i);
   if models{1}.fallback
       h3(i)=1; % Uniform fallback selected by the history-only preparation helper.
   else
       h3(i)=h3(i)/models{1}.Z;
   end
            
                m4=S4(i-1)-min(S4(max(i-5000,1):i-1));   
 
            % Select using only past log-capital drawups over W=5000.
            % If activation fails, the current factor stays at its initial value 1.
          [ii dd]=  max([m,m2,m3,m4,m_2,m_3,m_4]);
            if(ii>log(epsilon))
           diagnostics.selected(i)=dd;   
           switch dd
               case 1
             f2(i)=g(i);    
             choose=[choose; 1];
               case 2
                 f2(i)=h(i);
                  choose=[ choose; 2];
               case 3
                f2(i)=h2(i);
                 choose=[ choose; 3];
               case 4
                f2(i)=h3(i);
                 choose=[ choose; 4];
             case 5
                 f2(i)=h_(i);
                  choose=[ choose; 5];
               case 6
                f2(i)=h_2(i);
                 choose=[ choose; 6];
               case 7
                f2(i)=h_3(i);
                 choose=[ choose; 7];
           
           
           end     
              
            end   
            
          
             
% Update every candidate after selection; kNN factors are normalized.
                S3(i)= S3(i-1)+log(g(i));% Uniform-candidate log capital.
                S2(i)= S2(i-1)+log(h(i));S_2(i)= S_2(i-1)+log(h_(i));
                 S5(i)= S5(i-1)+log(h2(i));S_5(i)= S_5(i-1)+log(h_2(i));
                  S4(i)= S4(i-1)+log(h3(i));S_4(i)= S_4(i-1)+log(h_3(i));
% Linear capital can overflow or underflow; log capital is updated separately.
                S(i)=S(i-1)*f2(i);
                diagnostics.candidateValues(i,:)=[g(i),h(i),h2(i),h3(i),h_(i),h_2(i),h_3(i)];
                diagnostics.candidateLog(i,:)=[S3(i),S2(i),S5(i),S4(i),S_2(i),S_5(i),S_4(i)];
                logResults(i)=logResults(i-1)+log(f2(i));
  
% This comparison has no action; callers determine alarm crossings.
                if(S(i)>significance_l)
                
            end
            end
            
            
            
            results=S;
   
   
   
    end  
        

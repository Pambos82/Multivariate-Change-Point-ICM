function Tables = summarize_results(All)
%SUMMARIZE_RESULTS Summaries without changing the global RNG.
%   Z columns correspond to k=[5 10 15]. Candidate selection codes follow
%   combination. Startup observations 1 and 2 are excluded from
%   density/selection diagnostics. Delay means follow the original engine
%   (computed only on repetitions with post-change detection).

Performance = table();
NormalizationConstants = table();
SelectionFrequencies = table();
PerRepetition = table();
NumericalIssues = table();
paperNames = {'protein','particle','sensorless','credit'};
paperDetection = [.84,1,.99,.99];
paperFA = [.01,.01,0,.01];
paperMean = [83.9,15.1,60.7,83.5];
candidateNames = {'cautious/no bet','uniform','knn15','knn10','knn5','histogram5','histogram10','histogram15'};

for j=1:numel(All)
    item=All(j);
    paperIdx=find(strcmp(item.Dataset,paperNames),1);
    if isempty(paperIdx)
        pd=NaN; pf=NaN; pm=NaN;
    else
        pd=paperDetection(paperIdx); pf=paperFA(paperIdx); pm=paperMean(paperIdx);
    end
    run=item.Normalized;
    reps=numel(run.point);
    row=run.Results;
    row=addvars(row,{item.Dataset},{'normalized'},reps, ...
        'Before',1,'NewVariableNames',{'Dataset','Implementation','Repetitions'});
    row.PaperDetectionRate=pd;
    row.PaperFalseAlarmRate=pf;
    row.PaperMeanDelay=pm;
    row.DetectionMinusPaper=row.DetectionRate-pd;
    row.FalseAlarmMinusPaper=row.FalseAlarmRate-pf;
    row.MeanDelayMinusPaper=row.MeanDelay-pm;
    Performance=[Performance;row]; %#ok<AGROW>

    row=table(repmat({item.Dataset},reps,1),repmat({'normalized'},reps,1),(1:reps)', ...
        run.point>0,run.FalseAlarm,run.point, ...
        'VariableNames',{'Dataset','Implementation','Repetition','Detected','FalseAlarm','DelayZeroIfMiss'});
    PerRepetition=[PerRepetition;row]; %#ok<AGROW>
    row=table({item.Dataset},{'normalized'},sum(isinf(run.Martingale(:))), ...
        sum(isnan(run.Martingale(:))),sum(run.Martingale(:)==0), ...
        sum(~isfinite(run.LogMartingale(:))), ...
        'VariableNames',{'Dataset','Implementation','InfiniteLinearValues','NaNLinearValues','ZeroLinearValues','NonfiniteLogValues'});
    NumericalIssues=[NumericalIssues;row]; %#ok<AGROW>

    T=size(run.Martingale,2);
    selected=NaN(reps,T);
    for r=1:reps
        selected(r,:)=run.Diagnostics{r}.selected(:)';
    end
    assert(all(ismember(selected(:),0:7)), 'Unexpected candidate selection code.');
    phases={'all','pre_change','post_change'};
    phaseIndices={3:T,3:item.Parameters.nu,(item.Parameters.nu+1):T};
    for phase=1:3
        decisions=selected(:,phaseIndices{phase});
        denominator=numel(decisions);
        active=sum(decisions(:)>0);
        for code=0:7
            count=sum(decisions(:)==code);
            activeFrequency=NaN;
            if code>0 && active>0, activeFrequency=count/active; end
            frequency=NaN;
            if denominator>0, frequency=count/denominator; end
            row=table({item.Dataset},{'normalized'},{phases{phase}},code,{candidateNames{code+1}}, ...
                count,denominator,active,frequency,activeFrequency, ...
                'VariableNames',{'Dataset','Implementation','Phase','CandidateCode','Candidate', ...
                'Count','DecisionCount','ActiveBetCount','FrequencyAllDecisions','FrequencyAmongActiveBets'});
            SelectionFrequencies=[SelectionFrequencies;row]; %#ok<AGROW>
        end
    end
    zAll=NaN(reps*(T-2),3);
    fallbackAll=false(size(zAll));
    for r=1:reps
        index=(r-1)*(T-2)+(1:(T-2));
        zAll(index,:)=run.Diagnostics{r}.Z(3:end,:);
        fallbackAll(index,:)=run.Diagnostics{r}.fallback(3:end,:);
    end
    ks=[5,10,15];
    for c=1:3
        z=zAll(:,c);
        f=z(isfinite(z));
        q=NaN(1,9); avg=NaN; sd=NaN;
        if ~isempty(f)
            q=quantile(f,[0 .01 .05 .25 .5 .75 .95 .99 1]);
            avg=mean(f); sd=std(f);
        end
        row=table({item.Dataset},ks(c),numel(z),numel(f),sum(isnan(z)),sum(isinf(z)), ...
            sum(fallbackAll(:,c)),q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),avg,sd, ...
            'VariableNames',{'Dataset','K','DecisionCount','FiniteCount','NaNCount','InfiniteCount', ...
            'FallbackCount','Min','Q01','Q05','Q25','Median','Q75','Q95','Q99','Max','Mean','Std'});
        NormalizationConstants=[NormalizationConstants;row]; %#ok<AGROW>
    end
end
Tables=struct('Performance',Performance,'NormalizationConstants',NormalizationConstants, ...
    'SelectionFrequencies',SelectionFrequencies,'PerRepetition',PerRepetition, ...
    'NumericalIssues',NumericalIssues);
end

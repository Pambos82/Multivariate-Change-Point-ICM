function [Performance, outputDir] = run_full_figure_sweeps(options)
%RUN_FULL_FIGURE_SWEEPS Resumable full Figure 1/2 experiment.
%   [P,D] = run_full_figure_sweeps('DataDir',dataDir,'Workers',16)
%   evaluates K=1:100 and five ellipsoid lambda settings, plus the distance
%   NCM once per K. The 2,400 distinct jobs use normalized betting, epsilon
%   100, and 100 repetitions. Figure 1 uses the lambda=1e-3 subset.
%
%   Pass the returned D as 'OutputRoot',D to resume. Every finished job is
%   checked against the immutable manifest before reuse. Source/input hashes,
%   explicit RNG algorithms and transformations, p-values, decisions and
%   diagnostic summaries are retained. Incomplete temporary files are ignored.
%   Original scripts and existing results are never overwritten.
%
%   To test or run a bounded batch, use DatasetNames/KList/LambdaList/Reps
%   and MaxJobs. MaxJobs only limits newly scheduled jobs in this invocation;
%   it is not a scientific parameter and can change when resuming. Workers
%   can also change. An existing parallel pool must have the requested size.
%
%   All four datasets are prepared in historical Figure 2 order, including
%   when only a subset is selected. The scalar covariance behavior at NCM
%   K=1 is deliberately preserved. No full-grid validity is implied until
%   progress.json reports Status='complete' and all checkpoints exist.

arguments
    options.DataDir (1,:) char = ''
    options.DatasetNames = {'protein','particle','sensorless','credit'}
    options.KList (1,:) double = 1:100
    options.LambdaList (1,:) double = [1e-5 1e-4 1e-3 1e-2 1e-1]
    options.IncludeDistance (1,1) logical = true
    options.Reps (1,1) double = 100
    options.Epsilon (1,1) double = 100
    options.Workers (1,1) double = 16
    options.OutputRoot (1,:) char = ''
    options.MaxJobs (1,1) double = Inf
end
validateattributes(options.KList,{'double'},{'nonempty','finite','integer','positive','<=',4096});
validateattributes(options.LambdaList,{'double'},{'nonempty','finite','positive'});
validateattributes(options.Reps,{'double'},{'finite','integer','positive'});
validateattributes(options.Epsilon,{'double'},{'finite','positive'});
validateattributes(options.Workers,{'double'},{'finite','integer','positive'});
assert(isinf(options.MaxJobs) && options.MaxJobs>0 || ...
    isfinite(options.MaxJobs) && options.MaxJobs>=0 && fix(options.MaxJobs)==options.MaxJobs, ...
    'MaxJobs must be a nonnegative integer or Inf.');
options.KList = unique(options.KList,'stable');
options.LambdaList = unique(options.LambdaList,'stable');
validNames = {'protein','particle','sensorless','credit'};
options.DatasetNames = cellstr(string(options.DatasetNames));
assert(~isempty(options.DatasetNames) && all(ismember(options.DatasetNames,validNames)), ...
    'DatasetNames must contain protein, particle, sensorless, or credit.');
options.DatasetNames = validNames(ismember(validNames,options.DatasetNames));
repoDir = fileparts(mfilename('fullpath'));
addpath(repoDir,'-begin');
if isempty(options.DataDir), options.DataDir = getenv('COPA_DATA_DIR'); end
if isempty(options.DataDir), options.DataDir = fullfile(repoDir,'datasets'); end
options.DataDir = canonical_path(options.DataDir);
callerStream = RandStream.getGlobalStream;
callerState = callerStream.State;
restoreRng = onCleanup(@() restore_stream(callerStream,callerState)); %#ok<NASGU>

scientific = rmfield(options,{'Workers','OutputRoot','MaxJobs'});
scientific.Seed = 42;
scientific.Threshold = 100;
scientific.TRSize = 4096;
scientific.Nu = 256;
scientific.SigmaDelta = 0.1;
scientific.LookBackWindow = 5000;
scientific.ExperimentGenerator = 'threefry';
scientific.ExperimentNormalTransform = 'Inversion';
scientific.WorkerMaxNumCompThreads = 1;
scientific.PreprocessingGenerator = 'mt19937ar';
scientific.PreprocessingNormalTransform = 'Ziggurat';
scientific.PreprocessingOrder = validNames;
scientific.K1CovarianceBehavior = 'Preserve original scalar cov(single row,1)';
scientific.MatlabVersion = version;
sourceNames = {'run_full_figure_sweeps.m', ...
    'run_knn_cm_experiment.m','combination.m', ...
    'knn_betting_prepare.m','knn_ellipsoid_pdf.m', ...
    'pvalues.m','summarize_results.m'};
inputNames = {'CASP.mat','MiniBooNEPID.mat','Sensorlessdrivediagnosis.mat','creditcard_X.mat'};
sourceManifest = hash_files(repoDir,sourceNames);
inputManifest = hash_files(options.DataDir,inputNames);
signature = hash_bytes(unicode2native(jsonencode(struct('Scientific',scientific, ...
    'Sources',sourceManifest,'Inputs',inputManifest)),'UTF-8'));
jobs = make_jobs(options,validNames);

if isempty(options.OutputRoot)
    parentDir = fullfile(repoDir,'results','figure_sweeps');
    if ~isfolder(parentDir), mkdir(parentDir); end
    stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
    outputDir = fullfile(parentDir,['full_' stamp]);
    assert(~isfolder(outputDir),'Timestamp collision: choose a different OutputRoot.');
else
    outputDir = canonical_path(options.OutputRoot);
end
metadataPath = fullfile(outputDir,'metadata.mat');
if isfile(metadataPath)
    saved = load(metadataPath,'Metadata','jobs','K1Diagnostic');
    Metadata = saved.Metadata;
    assert(strcmp(Metadata.Signature,signature), ...
        'Resume rejected: scientific options, MATLAB version, source code or input hashes changed.');
    assert(isequaln(saved.jobs,jobs),'Resume rejected: job plan differs.');
    K1Diagnostic = saved.K1Diagnostic;
    currentPrepared = hash_files(fullfile(outputDir,'prepared'),{Metadata.PreparedFiles.Name});
    assert(isequaln(currentPrepared,Metadata.PreparedFiles), ...
        'Resume rejected: prepared data changed.');
    fprintf('Validated immutable manifest for resume: %s\n',outputDir);
else
    if isfolder(outputDir)
        listing = dir(outputDir);
        assert(all(ismember({listing.name},{'.','..'})), ...
            'OutputRoot is nonempty without a manifest; choose a fresh directory.');
    else
        mkdir(outputDir);
    end
    mkdir(fullfile(outputDir,'jobs'));
    mkdir(fullfile(outputDir,'prepared'));
    mkdir(fullfile(outputDir,'errors'));
    stream = RandStream('mt19937ar','Seed',42,'NormalTransform','Ziggurat');
    RandStream.setGlobalStream(stream);
    datasets = load_sweep_datasets(options.DataDir);
    preprocessingEndState = rng;
    K1Diagnostic = table();
    preparedNames = cell(1,numel(datasets));
    for d=1:numel(datasets)
        X = datasets(d).X; %#ok<NASGU>
        preparedNames{d} = [datasets(d).Name '.mat'];
        save(fullfile(outputDir,'prepared',preparedNames{d}),'X','-v7.3');
        covariance = cov(X(1,:),1);
        K1Diagnostic = [K1Diagnostic; table({datasets(d).Name},size(X,2), ...
            isscalar(covariance),covariance,'VariableNames', ...
            {'Dataset','Dimension','CovarianceIsScalar','CovarianceValue'})]; %#ok<AGROW>
    end
    clear datasets X;
    Metadata = struct('SchemaVersion',1,'CreatedUTC',utc_now(), ...
        'Scientific',scientific,'Signature',signature,'Sources',sourceManifest, ...
        'Inputs',inputManifest,'PreparedFiles',hash_files(fullfile(outputDir,'prepared'),preparedNames), ...
        'PreprocessingEndState',preprocessingEndState,'TotalJobs',numel(jobs), ...
        'Figure1HistoricalJitterSeedKnown',false, ...
        'Scope','Full requested grid; completion must be checked in progress.json', ...
        'CheckpointContents','P-values, per-repetition decisions, results and diagnostic summaries');
    save(metadataPath,'Metadata','jobs','K1Diagnostic','-v7');
    write_json_atomic(fullfile(outputDir,'metadata.json'),Metadata);
    writetable(struct2table(jobs),fullfile(outputDir,'PlannedJobs.csv'));
    writetable(K1Diagnostic,fullfile(outputDir,'K1CovarianceDiagnostic.csv'));
end

completed = false(1,numel(jobs));
for j=1:numel(jobs)
    path = checkpoint_path(outputDir,jobs(j));
    if isfile(path)
        validate_checkpoint(path,jobs(j),Metadata);
        completed(j)=true;
    end
end
pending = find(~completed);
scheduled = pending(1:min(numel(pending),options.MaxJobs));
fprintf('Jobs: %d validated complete / %d planned; scheduling %d now.\n', ...
    sum(completed),numel(jobs),numel(scheduled));
runTimer = tic;
failed = [];
write_progress(outputDir,Metadata,completed,failed,0,numel(scheduled),toc(runTimer),'running');

if ~isempty(scheduled)
    pool = gcp('nocreate');
    if isempty(pool)
        cluster = parcluster('Processes');
        % Honor the requested pool size without changing the saved profile.
        if cluster.NumWorkers < options.Workers
            cluster.NumWorkers = options.Workers;
        end
        pool = parpool(cluster,options.Workers);
    else
        assert(pool.NumWorkers==options.Workers, ...
            'Existing pool has %d workers, but Workers=%d. Reuse that size or close the pool first.', ...
            pool.NumWorkers,options.Workers);
    end
    nSlots = min(pool.NumWorkers,numel(scheduled));
    futures = parallel.FevalFuture.empty(1,0);
    futureMap = containers.Map('KeyType','double','ValueType','any');
    cancelOutstanding = onCleanup(@() cancel_pending(futureMap)); %#ok<NASGU>
    futureJobs = zeros(1,nSlots);
    next=1;
    for slot=1:nSlots
        id=scheduled(next); next=next+1; futureJobs(slot)=id;
        futures(slot)=parfeval(pool,@execute_job,1,repoDir,outputDir,jobs(id),Metadata);
        futureMap(id)=futures(slot);
    end
    finishedThisInvocation=0;
    while ~isempty(futures)
        [slot,result]=fetchNext(futures);
        id=futureJobs(slot);
        remove(futureMap,id);
        assert(result.JobNumber==id,'Worker returned an unexpected job number.');
        finishedThisInvocation=finishedThisInvocation+1;
        if result.Success
            validate_checkpoint(checkpoint_path(outputDir,jobs(id)),jobs(id),Metadata);
            completed(id)=true;
            fprintf('[%s] Complete %d/%d: job %d %s K=%d %s lambda=%.0e (%.1f s)\n', ...
                utc_now(),sum(completed),numel(jobs),id,jobs(id).Dataset,jobs(id).K, ...
                jobs(id).NCM,jobs(id).Lambda,result.ElapsedSeconds);
        else
            failed(end+1)=id; %#ok<AGROW>
            fprintf(2,'Job %d failed; checkpoint not accepted: %s\n',id,result.Error);
            write_json_atomic(fullfile(outputDir,'errors',sprintf('job_%04d_%s.json',id, ...
                char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS')))),result);
        end
        write_progress(outputDir,Metadata,completed,failed,finishedThisInvocation, ...
            numel(scheduled),toc(runTimer),'running');
        if next<=numel(scheduled)
            id=scheduled(next); next=next+1; futureJobs(slot)=id;
            futures(slot)=parfeval(pool,@execute_job,1,repoDir,outputDir,jobs(id),Metadata);
            futureMap(id)=futures(slot);
        else
            futures(slot)=[];
            futureJobs(slot)=[];
        end
    end
end
[Performance,SelectionFrequencies,NormalizationConstants,NumericalIssues,JobIndex] = ...
    aggregate_checkpoints(outputDir,jobs,completed);
if all(completed), status='complete'; elseif isempty(failed), status='partial'; else, status='failed'; end
write_progress(outputDir,Metadata,completed,failed,numel(scheduled),numel(scheduled),toc(runTimer),status);
save(fullfile(outputDir,'summary.tmp.mat'),'Performance','SelectionFrequencies', ...
    'NormalizationConstants','NumericalIssues','JobIndex','Metadata','K1Diagnostic','-v7');
movefile(fullfile(outputDir,'summary.tmp.mat'),fullfile(outputDir,'summary.mat'),'f');
fprintf('Figure sweep invocation %s: %d/%d completed.\n%s\n',status,sum(completed),numel(jobs),outputDir);
assert(isempty(failed),'Some jobs failed. See errors/ and resume after resolving the cause.');
end

function result = execute_job(repoDir,outputDir,job,metadata)
timer = tic;
result = struct('JobNumber',job.Number,'Success',false,'ElapsedSeconds',NaN,'Error','');
try
    addpath(repoDir,'-begin');
    maxNumCompThreads(metadata.Scientific.WorkerMaxNumCompThreads);
    assert(isequaln(hash_files(repoDir,{metadata.Sources.Name}),metadata.Sources), ...
        'Scientific source changed during the run; refusing an untracked computation.');
    path=checkpoint_path(outputDir,job);
    if isfile(path)
        validate_checkpoint(path,job,metadata);
        result.Success=true; result.ElapsedSeconds=toc(timer); return;
    end
    X = prepared_dataset(fullfile(outputDir,'prepared'),job.Dataset,metadata.Signature);
    stream=RandStream('threefry','Seed',42,'NormalTransform','Inversion');
    RandStream.setGlobalStream(stream);
    Job=job;
    Job.RngBefore=rng;
    Job.NormalTransform=stream.NormalTransform;
    % Capture the engine's routine repetition output; the controller reports
    % one concise completion message per job. Calculations are unchanged.
    evalc('datasetOut = run_knn_cm_experiment(X, ''reps'', metadata.Scientific.Reps, ''seed'', 42, ''TR_size'', 4096, ''nu'', 256, ''K'', job.K, ''sigma_delta'', 0.1, ''epsilon'', job.Epsilon, ''significance_l'', 100, ''cov_reg'', job.Lambda, ''ncm_type'', job.NCM, ''implementation'', ''normalized'');');
    datasetOut.Dataset=job.Dataset;
    Job.RngAfter=rng;
    actualStream=RandStream.getGlobalStream;
    assert(strcmp(actualStream.Type,stream.Type) && strcmp(actualStream.NormalTransform,'Inversion'), ...
        'Experiment changed the explicit RNG configuration.');
    Tables=summarize_results(datasetOut);
    % Remove Table 1 references: these are sweep streams, not Table 1 runs.
    Tables=rmfield(Tables,'Performance');
    jobPerformance=addvars(datasetOut.Normalized.Results,{job.Dataset},{'normalized'}, ...
        job.K,job.Lambda,{job.NCM},job.Epsilon,'Before',1,'NewVariableNames', ...
        {'Dataset','Implementation','K','Lambda','NCM','Epsilon'});
    datasetOut.Normalized=rmfield(datasetOut.Normalized,{'Martingale','LogMartingale','Diagnostics'});
    Job.ElapsedSeconds=toc(timer);
    Completion=struct('SchemaVersion',1,'Signature',metadata.Signature, ...
        'CompletedUTC',utc_now(),'Repetitions',metadata.Scientific.Reps,'StreamLength',512);
    temporary=[path(1:end-4) '.inprogress.mat'];
    save(temporary,'datasetOut','Job','Tables','jobPerformance','Completion','-v7');
    validate_checkpoint(temporary,job,metadata);
    assert(~isfile(path),'Completed checkpoint appeared during execution; refusing to overwrite it.');
    movefile(temporary,path);
    result.Success=true;
catch failure
    result.Error=getReport(failure,'extended','hyperlinks','off');
end
result.ElapsedSeconds=toc(timer);
end

function X = prepared_dataset(directory,name,signature)
% Each process loads a dataset at most once per run. Four datasets together
% fit comfortably in memory; they are not retransmitted with every job.
persistent savedDirectory savedSignature cached
if isempty(savedDirectory) || ~strcmp(savedDirectory,directory) || ~strcmp(savedSignature,signature)
    savedDirectory=directory; savedSignature=signature; cached=struct();
end
if ~isfield(cached,name)
    item=load(fullfile(directory,[name '.mat']),'X');
    cached.(name)=item.X;
end
X=cached.(name);
end

function validate_checkpoint(path,expected,metadata)
item=load(path,'Job','Completion','datasetOut','jobPerformance','Tables');
assert(all(isfield(item,{'Job','Completion','datasetOut','jobPerformance','Tables'})), ...
    'Checkpoint is incomplete: %s',path);
assert(strcmp(item.Completion.Signature,metadata.Signature) && item.Completion.SchemaVersion==1, ...
    'Checkpoint signature differs: %s',path);
names=fieldnames(expected);
for n=1:numel(names)
    assert(isfield(item.Job,names{n}) && isequaln(item.Job.(names{n}),expected.(names{n})), ...
        'Checkpoint job definition differs: %s',path);
end
reps=metadata.Scientific.Reps;
assert(isequal(size(item.datasetOut.PValues),[reps 512]) && ...
    all(isfinite(item.datasetOut.PValues(:))) && ...
    all(item.datasetOut.PValues(:)>=0 & item.datasetOut.PValues(:)<=1), ...
    'Checkpoint has invalid p-value streams: %s',path);
run=item.datasetOut.Normalized;
assert(~isfield(item.datasetOut,'Original') && isequal(size(run.point),[reps 1]) && ...
    isequal(size(run.FalseAlarm),[reps 1]) && all(ismember(run.point,0:256)) && ...
    islogical(run.FalseAlarm) && height(run.Results)==1 && height(item.jobPerformance)==1, ...
    'Checkpoint has invalid normalized decisions: %s',path);
assert(isequaln(run.Results,item.jobPerformance(:,7:end)), ...
    'Checkpoint result tables disagree: %s',path);
assert(height(item.Tables.NormalizationConstants)==3 && height(item.Tables.SelectionFrequencies)==24, ...
    'Checkpoint diagnostics are incomplete: %s',path);
parameters=item.datasetOut.Parameters;
assert(parameters.reps==reps && parameters.K==expected.K && parameters.cov_reg==expected.Lambda && ...
    parameters.epsilon==expected.Epsilon && parameters.significance_l==100 && parameters.TR_size==4096 && ...
    parameters.nu==256 && parameters.sigma_delta==0.1 && strcmp(parameters.ncm_type,expected.NCM) && ...
    strcmp(parameters.implementation,'normalized') && strcmp(item.datasetOut.Dataset,expected.Dataset), ...
    'Checkpoint engine parameters differ: %s',path);
delays=run.point(run.point>0);
if isempty(delays), meanDelay=NaN; medianDelay=NaN; stdDelay=NaN;
else, meanDelay=mean(delays); medianDelay=median(delays); stdDelay=std(delays); end
expectedResults=[mean(run.point>0),mean(run.FalseAlarm),mean(run.point==0),medianDelay,meanDelay,stdDelay];
assert(isequaln(table2array(run.Results),expectedResults), ...
    'Checkpoint summaries do not agree with its saved decisions: %s',path);
assert(all(item.Tables.NormalizationConstants.DecisionCount==reps*510) && ...
    all(item.Tables.NormalizationConstants.FallbackCount>=0) && ...
    height(item.Tables.NumericalIssues)==1 && height(item.Tables.PerRepetition)==reps, ...
    'Checkpoint diagnostic counts are incomplete: %s',path);
end

function jobs = make_jobs(options,validNames)
jobs=struct('Number',{},'Dataset',{},'K',{},'Lambda',{},'NCM',{},'Epsilon',{});
for d=1:numel(validNames)
    if ~ismember(validNames{d},options.DatasetNames), continue; end
    for k=options.KList
        for l=1:(numel(options.LambdaList)+double(options.IncludeDistance))
            if l<=numel(options.LambdaList)
                ncm='ellipsoid_pdf'; lambda=options.LambdaList(l);
            else
                ncm='knn_distance'; lambda=1e-3;
            end
            jobs(end+1)=struct('Number',numel(jobs)+1,'Dataset',validNames{d}, ...
                'K',k,'Lambda',lambda,'NCM',ncm,'Epsilon',options.Epsilon); %#ok<AGROW>
        end
    end
end
end

function [P,S,Z,N,J] = aggregate_checkpoints(directory,jobs,completed)
P=table(); S=table(); Z=table(); N=table(); J=table();
for id=find(completed)
    item=load(checkpoint_path(directory,jobs(id)),'Job','Tables','jobPerformance');
    job=item.Job;
    P=[P;item.jobPerformance]; %#ok<AGROW>
    S=[S;tag_table(item.Tables.SelectionFrequencies,job)]; %#ok<AGROW>
    Z=[Z;tag_table(item.Tables.NormalizationConstants,job)]; %#ok<AGROW>
    N=[N;tag_table(item.Tables.NumericalIssues,job)]; %#ok<AGROW>
    J=[J;table(job.Number,{job.Dataset},job.K,job.Lambda,{job.NCM},job.Epsilon, ...
        job.ElapsedSeconds,{sprintf('jobs/job_%04d.mat',job.Number)},'VariableNames', ...
        {'Job','Dataset','K','Lambda','NCM','Epsilon','ElapsedSeconds','File'})]; %#ok<AGROW>
end
tables={P,S,Z,N,J};
names={'Performance','SelectionFrequencies','NormalizationConstants','NumericalIssues','JobIndex'};
for i=1:numel(tables)
    temporary=fullfile(directory,[names{i} '.tmp.csv']);
    writetable(tables{i},temporary);
    movefile(temporary,fullfile(directory,[names{i} '.csv']),'f');
end
end

function T=tag_table(T,job)
if isempty(T), return; end
n=height(T);
T=addvars(T,repmat(job.K,n,1),repmat(job.Lambda,n,1),repmat({job.NCM},n,1), ...
    repmat(job.Epsilon,n,1),'Before',1,'NewVariableNames',{'NCMK','Lambda','NCM','Epsilon'});
end

function write_progress(directory,metadata,completed,failed,finished,scheduled,elapsed,status)
remaining=numel(completed)-sum(completed);
estimatedRemainingSeconds=NaN;
if finished>0, estimatedRemainingSeconds=elapsed/finished*remaining; end
if remaining==0, estimatedRemainingSeconds=0; end
progress=struct('UpdatedUTC',utc_now(),'Status',status,'Signature',metadata.Signature, ...
    'CompletedJobs',sum(completed),'TotalJobs',numel(completed),'CompletedJobNumbers',find(completed), ...
    'FailedThisInvocation',failed,'FinishedThisInvocation',finished,'ScheduledThisInvocation',scheduled, ...
    'InvocationElapsedSeconds',elapsed,'EstimatedRemainingSeconds',estimatedRemainingSeconds, ...
    'EstimateNote','Observed wall time per completed job; includes pool startup and may vary by dataset.', ...
    'Complete',all(completed));
write_json_atomic(fullfile(directory,'progress.json'),progress);
end

function cancel_pending(futureMap)
pending=values(futureMap);
for i=1:numel(pending), cancel(pending{i}); end
end

function path=checkpoint_path(directory,job)
path=fullfile(directory,'jobs',sprintf('job_%04d.mat',job.Number));
end

function datasets=load_sweep_datasets(directory)
datasets=struct('Name',{'protein','particle','sensorless','credit'},'X',[]);
S=load(fullfile(directory,'CASP.mat'),'CASP'); datasets(1).X=S.CASP{:,:};
S=load(fullfile(directory,'MiniBooNEPID.mat'),'MiniBooNEPID');
X=S.MiniBooNEPID{:,:}; datasets(2).X=X+1e-3*randn(size(X));
S=load(fullfile(directory,'Sensorlessdrivediagnosis.mat'),'Sensorlessdrivediagnosis');
X=S.Sensorlessdrivediagnosis{:,:};
if size(X,2)==49, X=X(:,1:48); end
datasets(3).X=X+1e-3*randn(size(X));
S=load(fullfile(directory,'creditcard_X.mat'),'creditcard_X'); datasets(4).X=S.creditcard_X;
end

function manifest=hash_files(directory,names)
manifest=struct('Name',{},'SHA256',{});
for i=1:numel(names)
    path=fullfile(directory,names{i});
    assert(isfile(path),'Required file is missing: %s',path);
    manifest(i)=struct('Name',names{i},'SHA256',hash_file(path)); %#ok<AGROW>
end
end

function hash=hash_file(path)
digest=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0,'Cannot read %s',path);
closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    chunk=fread(fid,4*1024*1024,'*uint8');
    if ~isempty(chunk), digest.update(typecast(chunk(:),'int8')); end
end
hash=digest_hex(digest);
end

function hash=hash_bytes(bytes)
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(uint8(bytes(:)),'int8'));
hash=digest_hex(digest);
end

function hash=digest_hex(digest)
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function path=canonical_path(path)
path=char(java.io.File(path).getCanonicalPath());
end

function stamp=utc_now()
stamp=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
end

function restore_stream(stream,state)
RandStream.setGlobalStream(stream); stream.State=state;
end

function write_json_atomic(path,value)
temporary=[path '.tmp'];
fid=fopen(temporary,'w'); assert(fid>=0,'Cannot create %s',temporary);
try
    fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
    fclose(fid);
catch failure
    fclose(fid); rethrow(failure);
end
movefile(temporary,path,'f');
end

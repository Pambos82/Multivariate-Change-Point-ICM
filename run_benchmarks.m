function varargout = run_benchmarks(options)
%RUN_BENCHMARKS Four benchmark datasets with normalized kNN betting only.
%   [All,Tables,D] = run_benchmarks() runs 100 repetitions per dataset,
%   cautious epsilon=100 and alarm threshold=100. Tables.Performance holds
%   detection/false-alarm rates and conditional delay summaries. Other tables
%   report normalization constants, selection frequencies and numerical issues.
%   Outputs are saved in a NEW results/benchmarks directory on every run.
%   For a smoke test use run_benchmarks('Reps',1). To replay saved p-values,
%   supply ReplayFile pointing to a previous benchmark.mat.
%   Console output is only the eight-column final performance summary.
%   Explicit outputs retain All, the complete Tables struct, and outputDir.

arguments
    options.DataDir (1,:) char = ''
    options.Reps (1,1) double = 100
    options.Epsilon (1,1) double = 100
    options.Threshold (1,1) double = 100
    options.OutputRoot (1,:) char = ''
    options.ReplayFile (1,:) char = ''
end
nargoutchk(0,3);
validateattributes(options.Reps, {'numeric'}, {'integer','positive'});
validateattributes(options.Epsilon, {'numeric'}, {'real','finite','positive'});
validateattributes(options.Threshold, {'numeric'}, {'real','finite','positive'});
repoDir = fileparts(mfilename('fullpath'));
addpath(repoDir);
if isempty(options.DataDir), options.DataDir = getenv('COPA_DATA_DIR'); end
if isempty(options.DataDir), options.DataDir = fullfile(repoDir,'datasets'); end
if isempty(options.OutputRoot), options.OutputRoot = fullfile(repoDir,'results','benchmarks'); end
if ~isfolder(options.OutputRoot), mkdir(options.OutputRoot); end
stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
baseDir = fullfile(options.OutputRoot, ['benchmark_' stamp]);
outputDir = baseDir;
counter = 0;
while isfolder(outputDir)
    counter = counter + 1;
    outputDir = sprintf('%s_%d',baseDir,counter);
end
mkdir(outputDir);
Metadata = struct('CreatedUTC',char(datetime('now','TimeZone','UTC')), ...
    'MatlabVersion',version,'Options',options,'Generator','twister', ...
    'LookBackWindow',5000,'KnnHistoryLength',100,'HistogramHistoryLength',200, ...
    'StartupExcludedFromDiagnostics',2,'BettingImplementation','normalized');
started = tic;

if ~isempty(options.ReplayFile)
    saved = load(options.ReplayFile,'All');
    assert(isfield(saved,'All'), 'Replay file must contain All from a benchmark run.');
    All = struct([]);
    for j = 1:numel(saved.All)
        params = saved.All(j).Parameters;
        params.reps = size(saved.All(j).PValues,1);
        params.epsilon = options.Epsilon;
        params.significance_l = options.Threshold;
        params.implementation = 'normalized';
        params.pvalue_streams = saved.All(j).PValues;
        args = namedargs2cell(params);
        datasetOut = run_knn_cm_experiment([],args{:});
        datasetOut.Dataset = saved.All(j).Dataset;
        if j==1, All=datasetOut; else, All(j)=datasetOut; end %#ok<AGROW>
        save(fullfile(outputDir,[datasetOut.Dataset '.mat']),'datasetOut','Metadata','-v7.3');
    end
else
    % Same order as run_benchmarks: starting RNG, protein experiment, particle
    % jitter/experiment, sensorless jitter/experiment, credit experiment.
    % No reporting operation consumes this global random-number stream.
    rng(42,'twister');
    names = {'protein','particle','sensorless','credit'};
    for j = 1:numel(names)
        switch names{j}
            case 'protein'
                S = load(fullfile(options.DataDir,'CASP.mat'),'CASP');
                X = S.CASP{:,1:end};
            case 'particle'
                S = load(fullfile(options.DataDir,'MiniBooNEPID.mat'),'MiniBooNEPID');
                X = S.MiniBooNEPID{:,1:end};
                X = X + 1e-3*randn(size(X));
            case 'sensorless'
                S = load(fullfile(options.DataDir,'Sensorlessdrivediagnosis.mat'),'Sensorlessdrivediagnosis');
                X = S.Sensorlessdrivediagnosis{:,1:end};
                if size(X,2)==49, X = X(:,1:48); end
                X = X + 1e-3*randn(size(X));
            case 'credit'
                S = load(fullfile(options.DataDir,'creditcard_X.mat'),'creditcard_X');
                X = S.creditcard_X;
        end
        clear S;
        datasetOut = run_knn_cm_experiment(X, ...
            'reps',options.Reps,'seed',42,'TR_size',4096,'nu',256, ...
            'K',50,'sigma_delta',0.1,'epsilon',options.Epsilon, ...
            'significance_l',options.Threshold,'implementation','normalized');
        datasetOut.Dataset = names{j};
        if j==1, All=datasetOut; else, All(j)=datasetOut; end %#ok<AGROW>
        save(fullfile(outputDir,[datasetOut.Dataset '.mat']),'datasetOut','Metadata','-v7.3');
    end
end

Tables = summarize_results(All);
Metadata.ElapsedSeconds = toc(started);
save(fullfile(outputDir,'benchmark.mat'),'All','Tables','Metadata','-v7.3');
tableNames = fieldnames(Tables);
for j=1:numel(tableNames)
    writetable(Tables.(tableNames{j}),fullfile(outputDir,[tableNames{j} '.csv']));
end
write_benchmark_report(fullfile(outputDir,'benchmark.md'),Tables,Metadata);
displayColumns = {'Dataset','Repetitions','DetectionRate','FalseAlarmRate', ...
    'MissRate','MedianDelay','MeanDelay','StdDelay'};
disp(Tables.Performance(:,displayColumns));
if nargout > 0
    varargout = {All, Tables, outputDir};
end
end

function write_benchmark_report(path,Tables,Metadata)
fid=fopen(path,'w'); assert(fid>=0,'Cannot write benchmark report.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# kNN betting benchmark\n\n');
fprintf(fid,'Generated UTC: %s. Executed in MATLAB.\n\n',Metadata.CreatedUTC);
fprintf(fid,'Cautious betting; epsilon=%g, alarm threshold=%g, W=5000. ', ...
    Metadata.Options.Epsilon,Metadata.Options.Threshold);
fprintf(fid,'Reference size=4096, two segments of 256, NCM K=50, covariance diagonal adjustment=1e-3, seed=42.\n\n');
fprintf(fid,'| Dataset | Repetitions | Detection | False alarm | Mean delay | Median delay | SD delay |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|\n');
P=Tables.Performance;
for j=1:height(P)
    fprintf(fid,'| %s | %d | %.4f | %.4f | %.8g | %.8g | %.8g |\n', ...
        P.Dataset{j},P.Repetitions(j),P.DetectionRate(j),P.FalseAlarmRate(j), ...
        P.MeanDelay(j),P.MedianDelay(j),P.StdDelay(j));
end
fprintf(fid,'\nDelays are conditional on post-change detection. A false alarm is any strict threshold crossing in the first segment. ');
fprintf(fid,'The first two betting positions are startup. Runs with fewer than 100 repetitions are smoke tests or exploratory runs, not the full paper protocol.\n\n');
fprintf(fid,'## Published Table 1 reference\n\n');
fprintf(fid,'Performance.csv includes the published detection/false-alarm/mean-delay values and current-minus-paper differences. ');
fprintf(fid,'These references used the historical implementation without density scaling; previous audits reproduced the reported values at epsilon=10, whereas the manuscript states 100. ');
fprintf(fid,'This variant at epsilon=%g is not claimed to produce identical published numbers.\n\n',Metadata.Options.Epsilon);
fprintf(fid,'## Output files\n\n');
fprintf(fid,'- benchmark.mat: exact p-values, martingales, full diagnostics and settings.\n');
fprintf(fid,'- Performance.csv: detection, false alarms, delay statistics and published-reference differences.\n');
fprintf(fid,'- NormalizationConstants.csv: finite-Z distributions and fallback/nonfinite counts for k=5,10,15.\n');
fprintf(fid,'- SelectionFrequencies.csv: candidate frequencies overall, before and after the change.\n');
fprintf(fid,'- PerRepetition.csv: decisions, false alarms and delays (zero denotes no post-change detection).\n');
fprintf(fid,'- NumericalIssues.csv: linear overflows/underflows and nonfinite log counts.\n\n');
fprintf(fid,'Diagnostics exclude the two startup positions. A uniform fallback replaces a nonintegrable or numerically unusable density. ');
fprintf(fid,'Analytical antiderivatives are evaluated in double precision; numerical normalization alone does not establish an anytime-valid guarantee.\n');
end

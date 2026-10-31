function [Stream, Summary, outputDir] = run_real(options)
%RUN_REAL Recreate the Air Quality experiment with normalized betting only.
%   [STREAM,SUMMARY,D] = RUN_REAL('DataDir',localAirDir) reads independently
%   prepared airquality.mat and saves a new results/airquality run. The raw
%   Air Quality input is deliberately not included in this package; see
%   datasets/PREPROCESSING.md and PREPARE_AIRQUALITY_INPUT.
%
%   Preserves the historical timetable preprocessing, reference size 3000,
%   K=50, covariance diagonal adjustment=1e-3, rng(42,'twister'), epsilon=100,
%   W=5000 and alarm threshold=100. Only normalized kNN betting is used.
%   Figure 3 uses accumulated log values over the complete monitoring
%   period, including where the linear martingale overflows.
arguments
    options.DataDir (1,:) char = ''
    options.OutputRoot (1,:) char = ''
end
repoDir = fileparts(mfilename('fullpath'));
addpath(repoDir,'-begin');
if isempty(options.DataDir), options.DataDir = getenv('COPA_DATA_DIR'); end
if isempty(options.DataDir), options.DataDir = fullfile(repoDir,'datasets'); end
inputPath = fullfile(options.DataDir,'airquality.mat');
assert(isfile(inputPath),'COPA:MissingAirQuality', ...
    ['Air Quality is not distributed. Download the original XLSX for your ' ...
     'authorized use, call prepare_airquality_input, and pass its folder as DataDir.']);
if isempty(options.OutputRoot)
    options.OutputRoot = fullfile(repoDir,'results','airquality');
end

loaded = load(inputPath,'AirQualityUCI');
assert(isfield(loaded,'AirQualityUCI') && istimetable(loaded.AirQualityUCI), ...
    'COPA:AirQualitySchema','Expected an AirQualityUCI timetable.');
T = loaded.AirQualityUCI;
assert(width(T)>=14,'COPA:AirQualitySchema','Expected at least 14 timetable variables.');
% These positional operations intentionally match the historical driver.
T = T(:,2:14);
T(:,3) = [];
T{:,:}(T{:,:} == -200) = NaN;
T = rmmissing(T);
X = T{:,:};
rowTimes = T.Properties.RowTimes;
assert(isnumeric(X) && isreal(X) && all(isfinite(X),'all') && ...
    size(X,1)>3000 && size(X,2)==12,'COPA:AirQualityValues', ...
    'Expected more than 3000 complete observations with 12 numeric measurements.');
assert(isdatetime(rowTimes) && ~any(isnat(rowTimes)) && issorted(rowTimes), ...
    'COPA:AirQualityTimes','Expected ordered historical datetime RowTimes.');
assert(isequal(rowTimes,dateshift(rowTimes,'start','day')), ...
    'COPA:AirQualityTimes', ...
    'Preserve the historical date-only RowTimes; use prepare_airquality_input.');
week = floor(days(rowTimes-rowTimes(1))/7)+1;

if ~isfolder(options.OutputRoot), mkdir(options.OutputRoot); end
stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
baseDir = fullfile(options.OutputRoot,['run_' stamp]);
outputDir = baseDir;
counter = 0;
while isfolder(outputDir)
    counter = counter+1;
    outputDir = sprintf('%s_%d',baseDir,counter);
end
mkdir(outputDir);
Metadata = struct('CreatedUTC',char(datetime('now','TimeZone','UTC')), ...
    'MatlabVersion',version,'Implementation','normalized','Seed',42, ...
    'Generator','twister','ReferenceSize',3000,'K',50, ...
    'CovarianceRegularization',1e-3,'Epsilon',100,'W',5000,'Threshold',100, ...
    'KnnHistoryLength',100,'HistogramHistoryLength',200, ...
    'CleanRows',size(X,1),'Features',size(X,2), ...
    'FirstCleanTimestamp',char(rowTimes(1)),'LastCleanTimestamp',char(rowTimes(end)), ...
    'InputSHA256',file_hash(inputPath),'GroundTruthChangePoint',NaN, ...
    'Week33IsVisualReferenceOnly',true);
sourceNames = {'run_real.m','combination.m','knn_betting_prepare.m', ...
    'knn_ellipsoid_pdf.m','pvalues.m'};
Metadata.Sources = repmat(struct('File','','SHA256',''),numel(sourceNames),1);
for j = 1:numel(sourceNames)
    Metadata.Sources(j) = struct('File',sourceNames{j}, ...
        'SHA256',file_hash(fullfile(repoDir,sourceNames{j})));
end
started = tic;
rng(42,'twister');
Stream = run_stream(X,week);
Stream.rowTimes = rowTimes(3001:end);
Metadata.ElapsedSeconds = toc(started);
Metadata.FinalRNG = rng;

idx = Stream.alarm_index;
if isempty(idx), idx = NaN; end
Summary = table({'normalized'},100,numel(Stream.pval),idx,Stream.alarm_time, ...
    Stream.alarm_week,Stream.logM(end),max(Stream.logM),sum(isinf(Stream.M)), ...
    sum(isnan(Stream.M)),all(isfinite(Stream.logM)), ...
    min(Stream.week_test),max(Stream.week_test), ...
    'VariableNames',{'Implementation','Epsilon','MonitoredSamples', ...
    'FirstAlarmMonitoredIndex','FirstAlarmCleanDataIndex','FirstAlarmWeek', ...
    'FinalLogM','MaxLogM','LinearOverflowCount','LinearNaNCount','FiniteLogM', ...
    'MonitoringStartWeek','MonitoringEndWeek'});
[SelectionFrequencies,NormalizationConstants] = diagnostic_tables(Stream.diagnostics);
save(fullfile(outputDir,'airquality_results.mat'),'Stream','Summary', ...
    'Metadata','SelectionFrequencies','NormalizationConstants','-v7.3');
writetable(Summary,fullfile(outputDir,'Summary.csv'));
writetable(SelectionFrequencies,fullfile(outputDir,'SelectionFrequencies.csv'));
writetable(NormalizationConstants,fullfile(outputDir,'NormalizationConstants.csv'));
write_text(fullfile(outputDir,'provenance.json'),jsonencode(Metadata,'PrettyPrint',true));
plot_stream(Stream,outputDir);
report = sprintf(['# Air Quality experiment\n\n' ...
    'MATLAB. Seed 42 (Twister), reference 3000, K=50, covariance diagonal adjustment 0.001, ' ...
    'epsilon=100, W=5000, alarm threshold=100.\n\n' ...
    'The %d complete observations have 12 measurements; %d are monitored. ' ...
    'Figure 3 contains the full monitoring interval and uses independently accumulated log martingale values. ' ...
    'Linear overflow count: %d. First alarm: monitored index %g, clean-data index %g, week %g.\n\n' ...
    'The alarm uses the strict rule log(M)>log(100). Week 33 is a visual reference, ' ...
    'not an assumed ground-truth change point. No detection delay or false-alarm rate is inferred ' ...
    'from this single unlabeled stream. The result uses kNN betting densities scaled by their integrals and is a methodological variant, ' ...
    'not an exact recovery of the published curve.\n\n' ...
    'Normalization uses analytical antiderivatives in floating-point arithmetic. Uniform fallbacks ' ...
    'replace nonintegrable or numerically unusable histories. Numerical normalization alone does ' ...
    'not establish the assumptions needed for a complete anytime-valid guarantee.\n'], ...
    Metadata.CleanRows,numel(Stream.pval),sum(isinf(Stream.M)),idx, ...
    Stream.alarm_time,Stream.alarm_week);
write_text(fullfile(outputDir,'airquality_results.md'),report);
disp(Summary);
fprintf('Air Quality outputs saved to:\n%s\n',outputDir);
end

function out = run_stream(X,week)
% Fixed settings and original arithmetic, apart from normalized betting.
[~,d] = size(X);
TR = X(1:3000,:);
TEST = X(3001:end,:);
mu = mean(TR,1);
sigma = std(TR,0,1);
sigma(sigma<1e-12) = 1;
TRz = (TR-mu)./sigma;
TESTz = (TEST-mu)./sigma;
Mdl = createns(TRz,'NSMethod','kdtree');
Vd = pi^(d/2)/gamma(d/2+1);
ncms = zeros(size(TESTz,1),1);
for i = 1:size(TESTz,1)
    x = TESTz(i,:);
    idx = knnsearch(Mdl,x,'K',50);
    density = knn_ellipsoid_pdf(TRz,x,idx,Vd,3000,1e-3);
    ncms(i) = -log(max(density,1e-300));
end
pv = pvalues(ncms);
assert(all(isfinite(ncms)) && all(isfinite(pv) & pv>=0 & pv<=1), ...
    'COPA:AirQualityStream','Invalid NCMs or conformal p-values.');
[M,logM,diagnostics] = combination(pv',100,100);
assert(all(isfinite(logM)) && ~any(isnan(M)), ...
    'COPA:AirQualityMartingale','Unexpected nonfinite log values or NaN martingale.');
out = struct('M',M(:),'logM',logM(:),'pval',pv,'ncms',ncms, ...
    'week_test',week(3001:end),'diagnostics',diagnostics);
out.alarm_index = find(out.logM>log(100),1,'first');
if isempty(out.alarm_index)
    out.alarm_time = NaN;
    out.alarm_week = NaN;
else
    out.alarm_time = 3000+out.alarm_index;
    out.alarm_week = out.week_test(out.alarm_index);
end
end

function [selection,normalization] = diagnostic_tables(D)
selection = table();
names = [{'cautious/no bet'},D.candidateNames];
count = numel(D.selected)-2;
for code = 0:7
    n = sum(D.selected(3:end)==code);
    row = table(code,names(code+1),n,count,n/count, ...
        'VariableNames',{'CandidateCode','Candidate','Count','DecisionCount','Frequency'});
    selection = [selection;row]; %#ok<AGROW>
end
normalization = table();
for j = 1:3
    z = D.Z(3:end,j);
    finite = z(isfinite(z));
    q = nan(1,7); avg = NaN; sd = NaN;
    if ~isempty(finite)
        q = prctile(finite,[0 1 5 50 95 99 100]);
        avg = mean(finite); sd = std(finite);
    end
    row = table(D.k(j),numel(z),numel(finite),sum(D.fallback(3:end,j)), ...
        sum(isinf(z)),sum(isnan(z)),avg,sd,q(1),q(2),q(3),q(4),q(5),q(6),q(7), ...
        'VariableNames',{'K','DecisionCount','FiniteCount','FallbackCount','InfiniteCount', ...
        'NaNCount','Mean','Std','Min','Q01','Q05','Median','Q95','Q99','Max'});
    normalization = [normalization;row]; %#ok<AGROW>
end
end

function plot_stream(Stream,outputDir)
f = figure('Visible','off','Color','w','Position',[100 100 1100 600]);
cleanup = onCleanup(@() close(f)); %#ok<NASGU>
ax = axes(f);
plot(ax,Stream.week_test,Stream.logM,'LineWidth',1.5,'Color',[0.12 0.35 0.65], ...
    'DisplayName','Cautious martingale');
hold(ax,'on');
yline(ax,log(100),'k--','Alarm threshold: log(100)','HandleVisibility','off');
xline(ax,33,':','Week 33 reference','Color',[0.5 0.5 0.5], ...
    'HandleVisibility','off','LabelVerticalAlignment','top');
if ~isempty(Stream.alarm_index)
    idx = Stream.alarm_index;
    scatter(ax,Stream.week_test(idx),Stream.logM(idx),40,[0.78 0.24 0.12], ...
        'filled','DisplayName','First threshold crossing');
end
xlim(ax,[min(Stream.week_test),max(Stream.week_test)]);
xlabel(ax,'Week since first cleaned observation');
ylabel(ax,'Natural log martingale');
title(ax,'Air Quality: cautious betting, epsilon = 100');
grid(ax,'on'); box(ax,'on'); legend(ax,'Location','northwest');
exportgraphics(f,fullfile(outputDir,'figure3.png'),'Resolution',200);
savefig(f,fullfile(outputDir,'figure3.fig'));
exportgraphics(f,fullfile(outputDir,'figure3.pdf'),'ContentType','vector');
end

function hash = file_hash(path)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(path,'rb'); assert(fid>=0,'Cannot read file: %s',path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    bytes = fread(fid,1024*1024,'*uint8');
    if ~isempty(bytes), md.update(typecast(bytes(:),'int8')); end
end
hash = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function write_text(path,text)
fid = fopen(path,'w'); assert(fid>=0,'Cannot write file: %s',path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',text);
end

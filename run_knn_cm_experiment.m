function out = run_knn_cm_experiment(X, params)
%RUN_KNN_CM_EXPERIMENT Experiment engine using normalized cautious betting only.
%   Data sampling, standardization, NCMs and randomized p-values preserve the
%   research engine. Set pvalue_streams to replay a repetitions-by-time
%   matrix without resampling data. implementation is fixed to normalized.
%   The original unnormalized betting implementation is not required.

arguments
    X double
    params.reps (1,1) double = 100
    params.seed (1,1) double = 42
    params.TR_size (1,1) double = 4096
    params.nu (1,1) double = 256
    params.K (1,1) double = 50
    params.sigma_delta (1,1) double = 0.1
    params.epsilon (1,1) double = 100
    params.significance_l (1,1) double = 100
    params.standardize (1,1) logical = true
    params.change_type (1,:) char = 'mean_shift'
    params.cov_reg (1,1) double = 1e-3
    params.ncm_type (1,:) char = 'ellipsoid_pdf'
    params.implementation (1,:) char = 'normalized'
    params.pvalue_streams double = []
end

validateattributes(params.reps, {'numeric'}, {'integer','positive'});
validateattributes(params.nu, {'numeric'}, {'integer','positive'});
validateattributes(params.epsilon, {'numeric'}, {'finite','positive'});
validateattributes(params.significance_l, {'numeric'}, {'finite','positive'});
assert(strcmpi(params.implementation,'normalized'), ...
    'COPA:NormalizedOnly','This package supports normalized betting only.');
params.implementation = 'normalized';
replay = ~isempty(params.pvalue_streams);
T = 2 * params.nu;
if replay
    validateattributes(params.pvalue_streams, {'numeric'}, ...
        {'real','finite','>=',0,'<=',1,'size',[params.reps T]});
    streams = params.pvalue_streams;
else
    [N, d] = size(X);
    if N < params.TR_size + T
        error('Not enough rows: need %d, got %d.', params.TR_size + T, N);
    end
    % Deliberately matches the original engine, including its generator
    % inheritance. The benchmark runner sets rng(42,'twister') beforehand.
    rng(params.seed);
    Vd = pi^(d/2) / gamma(d/2 + 1);
    nRef = params.TR_size;
    lam = params.cov_reg;
    streams = NaN(params.reps,T);
end

out = struct();
out.Dataset = '';
out.Parameters = rmfield(params, 'pvalue_streams');
out.Replayed = replay;
out.Normalized = initialize_run(params.reps,T);

for rep = 1:params.reps
    if ~replay
        idx = randperm(N);
        TR = X(idx(1:params.TR_size), :);
        POOL = X(idx(params.TR_size+1:end), :);
        if params.standardize
            mu = mean(TR, 1);
            sigma = std(TR, 0, 1);
            sigma(sigma < 1e-12) = 1;
            TRz = (TR - mu) ./ sigma;
            POOLz = (POOL - mu) ./ sigma;
        else
            TRz = TR;
            POOLz = POOL;
        end
        Mdl = createns(TRz, 'NSMethod', 'kdtree');
        permY = randperm(size(POOLz,1));
        Y0 = POOLz(permY(1:params.nu), :);
        Y1 = POOLz(permY(params.nu+1:2*params.nu), :);
        if strcmpi(params.change_type, 'mean_shift')
            Delta = randn(1,d) * params.sigma_delta;
            Y1s = Y1 + Delta;
        elseif strcmpi(params.change_type, 'none')
            Y1s = Y1;
        else
            error('Unknown change_type: %s', params.change_type);
        end
        switch lower(params.ncm_type)
            case 'ellipsoid_pdf'
                idx0 = knnsearch(Mdl, Y0, 'K', params.K);
                idx1 = knnsearch(Mdl, Y1s, 'K', params.K);
                p0 = knn_ellipsoid_pdf(TRz, Y0, idx0, Vd, nRef, lam);
                p1 = knn_ellipsoid_pdf(TRz, Y1s, idx1, Vd, nRef, lam);
                ncms0 = -log(max(p0, 1e-300));
                ncms1 = -log(max(p1, 1e-300));
            case 'knn_distance'
                [~, dist0] = knnsearch(Mdl, Y0, 'K', params.K);
                [~, dist1] = knnsearch(Mdl, Y1s, 'K', params.K);
                ncms0 = dist0(:, end);
                ncms1 = dist1(:, end);
            otherwise
                error('Unknown ncm_type: %s', params.ncm_type);
        end
        pv = pvalues([ncms0; ncms1]);
        streams(rep,:) = pv(:)';
    end
    pv = streams(rep,:);
    % Betting must not consume the experiment random-number stream.
    stateBefore = rng;
    [M, logM, diagnostic] = combination(pv, params.epsilon, params.significance_l);
    out.Normalized = record_run(out.Normalized, rep, M, logM, diagnostic, params);
    assert(isequaln(stateBefore,rng), 'Betting changed the experiment RNG.');
end
out.PValues = streams;
out.Normalized.Results = summarize_run(out.Normalized);
end

function run = initialize_run(reps,T)
run = struct('point',zeros(reps,1), 'FalseAlarm',false(reps,1), ...
    'Martingale',NaN(reps,T), 'LogMartingale',NaN(reps,T), ...
    'Diagnostics',{cell(reps,1)}, 'Results',table());
end

function run = record_run(run,rep,M,logM,diagnostic,params)
M = M(:)'; logM = logM(:)';
assert(numel(M)==2*params.nu && numel(logM)==numel(M), 'Unexpected martingale length.');
assert(~any(isnan(M)), 'Martingale contains NaN.');
idxHit = find(M(params.nu+1:end) > params.significance_l, 1, 'first');
if ~isempty(idxHit), run.point(rep) = idxHit; end
run.FalseAlarm(rep) = any(M(1:params.nu) > params.significance_l);
run.Martingale(rep,:) = M;
run.LogMartingale(rep,:) = logM;
run.Diagnostics{rep} = diagnostic;
end

function result = summarize_run(run)
del = run.point(run.point > 0);
if isempty(del)
    medianDelay = NaN; meanDelay = NaN; sdDelay = NaN;
else
    medianDelay = median(del); meanDelay = mean(del); sdDelay = std(del);
end
result = table(mean(run.point > 0), mean(run.FalseAlarm), mean(run.point == 0), ...
    medianDelay, meanDelay, sdDelay, 'VariableNames', ...
    {'DetectionRate','FalseAlarmRate','MissRate','MedianDelay','MeanDelay','StdDelay'});
end

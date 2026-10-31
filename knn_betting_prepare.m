function model = knn_betting_prepare(history, k)
%KNN_BETTING_PREPARE Integrate the original one-dimensional kNN density.
%   MODEL = KNN_BETTING_PREPARE(HISTORY,K) depends only on past p-values.
%   For N >= 2 and m=min(N,K), the original density, away from zero
%   distances, is c/r_m(u), with c=min(N-1,K-1)/(2*N). Its integral over
%   [0,1] is computed by analytical logarithmic antiderivatives, in
%   floating-point arithmetic. There is no boundary reflection or clipping
%   of the radius, and no grid quadrature in this function.
%
%   If x is sorted, its m nearest neighbours form a contiguous block.
%   Block j is active between (x(j-1)+x(j+m-1))/2 and
%   (x(j)+x(j+m))/2, with outer endpoints 0 and 1. Within that block,
%   r_m(u)=max(u-x(j), x(j+m-1)-u). Splitting at the midpoint of the
%   two block endpoints gives the two logarithmic integrals below.
%
%   An atom repeated at least m times makes the integral infinite:
%   changing the value only at r_m(u)=0 (the original dist=1 rule) cannot
%   remove the adjacent 1/abs(u-x) singularity. Such histories use a
%   documented uniform betting fallback, NOT a normalized original kNN
%   density. The same fallback is used for N<2 or numerical failure.
%
%   MODEL.Z is the raw density's integral (Inf for divergence, 0 for
%   N=1, and NaN when N=0). MODEL.status and MODEL.fallback identify the
%   outcome. COMBINATION uses these prepared normalization constants.

validateattributes(history, {'numeric'}, {'real'}, mfilename, 'history');
if ~isvector(history) && ~isempty(history)
    error('knn_betting_prepare:InvalidHistory','History must be a vector.');
end
validateattributes(k, {'numeric'}, {'real','scalar','finite','integer','>=',2}, mfilename, 'k');
x = sort(double(history(:)));
if any(~isfinite(x) | x < 0 | x > 1)
    error('knn_betting_prepare:InvalidHistory', ...
        'Historical p-values must be finite and lie in [0,1].');
end
N = numel(x);
m = min(N,k);
model = struct('historySorted',x,'N',N,'k',double(k),'m',m, ...
    'coefficient',NaN,'Z',NaN,'status','uniform_small_sample', ...
    'fallback',true,'transitions',zeros(0,1));
if N == 0
    return;
end
model.coefficient = min(N-1,k-1)/(2*N);
if N == 1
    model.Z = 0;
    return;
end

% Block boundaries are also sufficient to evaluate the kth distance.
model.transitions = x(1:N-m) + (x(1+m:N)-x(1:N-m))/2;
first = x(1:N-m+1);
last = x(m:N);
if any(first == last)
    model.Z = Inf;
    model.status = 'uniform_divergent';
    return;
end
% Decide numerical fallbacks from history alone. The smallest possible
% mth-neighbour radius is half the narrowest m-point block span. If the
% corresponding raw peak cannot be represented, evaluating only the current
% point and substituting a bet afterward would change its integral.
minimumRadius = min(last-first)/2;
maximumRaw = model.coefficient/minimumRadius;
if minimumRadius <= 0 || ~isfinite(maximumRaw)
    model.status = 'uniform_numerical_failure';
    return;
end
left = [0; model.transitions];
right = [model.transitions; 1];
middle = first + (last-first)/2;

% Integral on [left,min(right,middle)] is log((last-left)/(last-end)).
% log1p retains precision for a narrow interval with a large radius.
lowEnd = min(right,middle);
lowWidth = max(0,lowEnd-left);
lowRadius = last-lowEnd;
highStart = max(left,middle);
highWidth = max(0,right-highStart);
highRadius = highStart-first;
lowActive = lowWidth > 0;
highActive = highWidth > 0;
terms = [log1p(lowWidth(lowActive)./lowRadius(lowActive)); ...
         log1p(highWidth(highActive)./highRadius(highActive))];
Z = model.coefficient * sum(terms);
model.Z = Z;
if ~isfinite(Z) || Z <= 0 || any(~isfinite(terms)) || ~isfinite(maximumRaw/Z)
    model.status = 'uniform_numerical_failure';
    return;
end
model.fallback = false;
model.status = 'analytic';
end

function p = knn_ellipsoid_pdf(TR, Q, idx, Vd, nRef, lam)
%TR  : nRef x d reference (standardized)
%Q   : m x d queries
%idx : m x K neighbor indices into TR

    [m, d] = size(Q);
    K = size(idx,2);
    p = zeros(m,1);

    for i = 1:m
        Xi = TR(idx(i,:), :);         % K x d neighbor points
        xq = Q(i,:);                  % 1 x d query

        % Local covariance of neighbors with lam added to the diagonal.
        % (center by neighbor mean; you could also center by xq)
        S = cov(Xi, 1) + lam * eye(d);

        % Cholesky for inverse + determinant
        % S should be positive definite after adding lam*I.
        R = chol(S);

        % Mahalanobis distances of neighbors to query
        D = Xi - xq;                  % K x d
        Z = D / R;                    % solves Z*R = D  => Z = D*inv(R)
        m2 = sum(Z.^2, 2);            % squared Mahalanobis distances
        t  = sqrt(max(m2));           % radius to cover the K neighbors

        % Ellipsoid volume: Vd * t^d * sqrt(det(S))
        % det(S) = (prod(diag(R))^2) so sqrt(det(S)) = prod(diag(R))
        sqrtDetS = prod(diag(R));
        vol = Vd * (max(t, 1e-12)^d) * sqrtDetS;

        % kNN density
        p(i) = K / (nRef * vol);
    end
end

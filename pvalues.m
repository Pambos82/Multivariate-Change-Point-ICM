function pvalue = pvalues(x)

    if isnumeric(x)
        ncms = x(:);
        n = length(ncms);

        pvalue = NaN(n,1);
        if n == 0
            return;
        end

        pvalue(1) = rand;

        for i = 2:n
            pvalue(i) = ( ...
                sum(ncms(1:i-1) > ncms(i)) + ...
                rand * sum(ncms(1:i) == ncms(i)) ) / i;
        end

    elseif iscell(x)
        out = x;
        n = numel(out);

        pvalue = NaN(n,1);

        for k = 1:n
            if isempty(out{k}) || ~isfield(out{k}, 'ncm') || isempty(out{k}.ncm) || all(isnan(out{k}.ncm))
                continue;
            end

            alpha = out{k}.ncm(:);
            valid = isfinite(alpha);
            alpha = alpha(valid);

            m = numel(alpha);
            if m == 0
                continue;
            end

            a_new = alpha(end);

            pvalue(k) = ( ...
                sum(alpha(1:end-1) > a_new) + ...
                rand * sum(alpha == a_new) ) / m;
        end

    elseif isstruct(x)
        n = numel(x);
        pvalue = NaN(n,1);

        for k = 1:n
            if ~isfield(x(k), 'ncm') || isempty(x(k).ncm) || all(isnan(x(k).ncm))
                continue;
            end

            alpha = x(k).ncm(:);
            valid = isfinite(alpha);
            alpha = alpha(valid);

            m = numel(alpha);
            if m == 0
                continue;
            end

            a_new = alpha(end);

            pvalue(k) = ( ...
                sum(alpha(1:end-1) > a_new) + ...
                rand * sum(alpha == a_new) ) / m;
        end

    else
        error('Input must be either a numeric vector, a cell array, or a struct array.');
    end
end

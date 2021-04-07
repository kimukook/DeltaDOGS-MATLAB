function x = normalize_bounds(x0, lb, ub)
    % Normalize the points to be within [0, 1] in each dimension.
    n = size(lb, 1);
    m = size(x0, 2);
    x = zeros(n, m);
    for i = 1 : n
        x(i, :) = (x0(i, :) - lb(i)) / (ub(i) - lb(i));
    end
end
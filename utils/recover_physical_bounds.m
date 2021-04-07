function x = recover_physical_bounds(x0, lb, ub)
    % Recover the point in physical constraints for function
    % evaluation.
    n = size(lb, 1);
    m = size(x0, 2);
    x = zeros(n, m);
    for i = 1 : n
        x(i, :) = x0(i, :) *(ub(i) - lb(i)) + lb(i);
    end
end
function xU = generate_support_points(lb, ub, n, xE)
    % Generate the support points (corner of the domain) without
    % repeated points in evaluated data set.
    xU = generate_bounds(lb, ub, n);
    non_repeat_index_list = [];
    for i = 1 : size(xU, 2)
       [dis, ~, ~] = mindis(xU(:, i), xE);
       if dis > 1e-6
           non_repeat_index_list = [non_repeat_index_list, i];
       end
    end
    xU = xU(:, non_repeat_index_list);
end
function xE = random_initial(n, m, Nm)
    % If there is no evaluated points as the input, DeltaDOGS
    % randomly select 2*dim points as the initial.
    xU = generate_bounds(zeros(n, 1), ones(n, 1), n);
    xE = round(rand(n, 1) * Nm) / Nm;
    while size(xE, 2) < m
        temp = round(rand(n, 1) * Nm) / Nm;
        [dis1, ~, ~] = mindis(temp, [xU, xE]);
        if dis1 > 1e-6
            xE = [xE, temp];
        end
    end
end
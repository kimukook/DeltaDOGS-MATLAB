% Name         :  Natural Polyharmonic Spline (NPS) interpolation
% Functionality:  Create the NPS interpolation for input x and output y
% Properties   :  xi -> sites
%                 yi -> function values
%              :  w  -> coefficients for residuals at each site
%              :  v  -> linear coefficients
%
% Author       :  Muhan Zhao
% Date         :  May. 16, 2019


classdef NPSInterpolation
    properties
        xi
        yi
        w
        v
        n
    end
    methods
        function inter_par = interpolateparametarization(inter_par, xi, yi)
        inter_par.xi = xi;
        inter_par.yi = yi;
        inter_par.n = size(xi,1);
        % polyharmonic spline interpolation
        N = size(xi,2); A = zeros(N,N);
        for ii = 1 : 1 : N
            for jj = 1 : 1 : N
                A(ii,jj) = ((xi(:,ii) - xi(:,jj))' * (xi(:,ii) - xi(:,jj))) ^ (3 / 2);
            end
        end
        V = [ones(1,N); xi];
        A = [A V'; V zeros(inter_par.n + 1, inter_par.n + 1)];
        wv = A \ [yi.'; zeros(inter_par.n+1,1)]; % solve the associated linear system
        inter_par.w = wv(1:N); 
        inter_par.v= wv(N+1:N+inter_par.n+1); 
        end
        
        function y = interpolate_eval(inter_par, x)
            S = bsxfun(@minus, inter_par.xi, x);
            y = inter_par.v' * [1; x] + inter_par.w' * sqrt(diag(S' * S)) .^ 3;
        end
        
        function g = interpolate_grad(inter_par, x)
            N = size(inter_par.xi, 2);
            g = zeros(inter_par.n, 1);
            for i = 1 : N
                X = x - inter_par.xi(:, i);
                g = g + 3 * inter_par.w(i) * X * norm(X);
            end
            g = g + inter_par.v(2:end);
        end 
    end
end

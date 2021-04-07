% Name         :  Delaunay-based Derivative-Free Optimization via Global
%                 Surrogate (DeltaDOGS) Matlab Code
%
% Functionality:  Optimize(Minimize) an unknown function via minimizing
%                 the surrogate model built upon Delaunay triangulation
%
% Example      :  To use the following data structure, try the following 
%                 2D schwefel example:
%                 ========================================================
%                 clear, clc, close all
%                 addpath(genpath(pwd));
%                 n = 2;
%                 xE = [0.2500    ,0.7500    ,0.5000;
%                       0.2500   , 0.5000  ,  0.87500];
%                 bounds = [zeros(n, 1), ones(n, 1)];
%                 surrogate_ = 'constant'; 
%                 % or use 'adaptive' to switch to constant surrogate model
%                 % if a target value is given as y0;
%                 func_eval = @(x) -sum((2*x) .* sin(sqrt(500*abs(x))));
%                 xmin = 0.8419 * ones(n, 1);
%                 y0 = -1.6759 * n;
%                 K = 3;
%                 mesh_size = 8;
%                 num_mesh_refine = 7;
%                 max_iter = 35;
%                 verbose = 1;
% 
%                 ddogs = DeltaDOGS;
%                 ddogs.initial(n, bounds, surrogate_, func_eval, xE, xmin, y0, K, num_mesh_refine, mesh_size, max_iter, verbose);
%                 ddogs.DeltaDogsOptimize;
% ========================================================
% Author       :  Muhan Zhao
% Institute    :  Mechanical and Aerospace Engineering, UC San Diego
% Date         :  Mar. 12, 2021

classdef DeltaDOGS < handle
    % This is skeleton code for Delaunay Triangulation search
    % Notice that we have '< handle' to change private properties of 
    % DeltaDOGS at each sub-function call. Otherwise this is a value class 
    % and do not change the inherit properties.
    properties
        % sites for constructing Delaunay triangulation
        xE      % Evaluated points.
        yE      % Function values associated with 'xE'.
        yRange  % Range of the objective function values
        xU      % Support points
        yU      % Discrete search function value of the support points
        xall    % stack xU and xE
        
        % Objective function
        func_eval  % The function evaluator
        
        % Physical upper and lower bounds
        physical_ub
        physical_lb
        
        % Mesh size info
        mesh_size            % The current mesh size
        num_mesh_refinement  % The maximum times of mesh refinement
        
        % Iteration information
        num_iter             % The number of iteration
        max_iter             % maximum number of iterations
        
        % Upper and lower bound of the normalized box domain
        ub  % all to be 1
        lb  % all to be 0
        
        % Delaunay triangulation info
        tri  % The index of Delauany triangulation
        
        % General information
        n              % The size of the input parameter space
        surrogate_type % The type of the surrogate model, 'constant' or 'adaptive'
        
        inter_par  % stores the parameters for interpolation
        y0         % The target value
        K          % The tradeoff parameter for constant K search model
        
        % Stores the associated function values at the center of Delaunay 
        % Triangulation
        ind_min % the index of the minimum of the evaluated point
        
        % Define the linear constraints, Ax <= b. if A and b are None type, set them to be the box domain constraints.
        Ain
        bin
        
        % Iterative sampling at each iteration
        xm              % minimizer of search model
        ym              % surrogate model value of xm
        xm_eval         % point to evaluate (quantized xm)
        refine_trigger  % trigger of mesh refinement
        
        % plot parameters: 
        obj_lim % obj_lim: is the objective function values range
        verbose % display the optimization information
        
        % result of optimization
        x % Global optimizer
    end
    
    methods
        function initial(self, n, bounds, surrogate_type, func_eval, xE, y0, K, num_mesh_refine, ms, max_iter, verbose)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is the initialization of Delaunay triangulation
            % optimization, the input information is declared by users 
            % given as follows:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % n                      :  The dimension of input parameter
            % bounds                 :  The physical upper and lower bounds
            % surrogate_type         :  The surrogate model type
            % func_eval              :  The function evaluator
            % xE                     :  The initial points to start,
            %                           including corners of the box domain
            % y0                     :  The target value
            % K                      :  The tradeoff parameter of constantK
            % num_mesh_size(optional):  The maximum time of mesh refinement
            % ms(optional)           :  The initial mesh size
            % max_iter(optional)     :  The maximum number of iterations
            % verbose(optional)      :  The display of output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin < 12
                verbose = 0;
            end
            if nargin < 11 || isempty(max_iter)
                max_iter = Inf;
            end
            if nargin < 10
                ms = 8;  % If no input of mesh size, initialize it with 8
            end
            if nargin < 9
                % If no input of mesh refinement, initialize it as 8
                num_mesh_refine = 8; 
            end
            self.num_iter = 0;
            self.n = n;
            self.func_eval = func_eval;
            self.physical_lb = bounds(:, 1);
            self.physical_ub = bounds(:, 2);
            self.lb = zeros(n, 1);
            self.ub = ones(n, 1);
            self.yE = zeros(1, size(self.xE, 2));
            self.surrogate_type = surrogate_type;
            self.num_mesh_refinement = num_mesh_refine;
            self.max_iter = max_iter;
            self.mesh_size = ms;
            self.verbose = verbose;
            if strcmp(self.surrogate_type, 'constant')
                self.K = K;
                self.y0 = [];
            elseif strcmp(self.surrogate_type, 'adaptive')
                self.K = [];
                self.y0 = y0;
            else
                disp("No surrogate model defined. Please set the property 'surrogate_type' to be either 'constant' or 'adaptive'. ")
            end
            self.yE = zeros(1, size(self.xE, 2));
            
            % Take in or generate the initial evaluated points set
            if size(xE, 2) > 0 
                self.xE = normalize_bounds(xE, self.physical_lb, self.physical_ub);
            else
                self.xE = random_initial(self.n, 2 * self.n, self.mesh_size);
            end
            
            % Generate the set of support points -> corner of the box
            % constraint.
            self.xU = generate_support_points(self.lb, self.ub, self.n, self.xE);
            
            % Define the linear constraints
            self.Ain = [eye(self.n); -eye(self.n)]; 
            self.bin = [ones(self.n, 1); zeros(self.n, 1)];
            
            % Evaluate the initial points
            for i = 1 : size(self.xE, 2)
                point_to_eval = recover_physical_bounds(self.xE(:, i), self.physical_lb, self.physical_ub);
                self.yE(i) = func_eval(point_to_eval);
            end
            
            % create the folder that stores all the figures
            if ~exist('figures', 'dir')
                mkdir('figures')
            end
        end
        
        function DeltaDogsOptimize(self)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is the main function which is called to optimize the
            % objective
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for kk = 1 : self.num_mesh_refinement
                for k = 1 : 100
                    % update the iteration number
                    self.num_iter = self.num_iter + 1;
                    if self.num_iter > self.max_iter
                        break
                    end
                    % reset the refine trigger
                    self.refine_trigger = 0;
                    self.DelaunaySearch;
                    if self.refine_trigger == 1
                        break
                    end
                end
                if self.num_iter > self.max_iter
                    break
                end
            end
            [~, ind] = min(self.yE);
            self.x = recover_physical_bounds(self.xE(:, ind), self.physical_lb, self.physical_ub);
        end
        
        function DelaunaySearch(self)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Iteratively build up the Delaunay triangulation and minimize
            % the surrogate model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create the interpolation
            interpolation = NPSInterpolation;
            self.inter_par = interpolation.interpolateparametarization(self.xE, self.yE);
            [~, self.ind_min] = min(self.yE);
            self.yRange = range(self.yE);
            
            % Construct the Delaunay triangulation and minimize the 
            % surrogate model
            if strcmp(self.surrogate_type, 'constant')
                constantDelaunayMinmization(self);
            elseif strcmp(self.surrogate_type, 'adaptive')
                adaptiveDelaunayMinimization(self);
            end
                
            % Generate the plot
            % plot_maker2D(self);
            % Update the evaluated data and proceeding
            DelaunayUpdate(self);
            if self.verbose == 1
                IterOptInfoOutput(self);
            end
        end
        
        function constantDelaunayMinmization(self)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the next point to evaluate through minimizing the 
            % constant surrogate model built upon Delaunay triangulation,
            % including evaluate the support points on the constraints.
            
            % Constant surrogate model: s(x) = p(x) - K * e(x)
            % Input :  Evaluated data points, xE and yE
            % Return:  Next point to evaluate, xm_eval
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compute the discrete search model values at the support point
            self.yU = zeros(size(self.xU, 2), 1);
            for ii = 1 : size(self.xU, 2)
                self.yU(ii) = self.inter_par.interpolate_eval(self.xU(:, ii)) - self.K * self.yRange * mindis(self.xU(:, ii), self.xE);
            end
            if size(self.xU, 2) > 0 && min(self.yU) < min(self.yE)
                [~, ind] = min(self.yU);
                self.xm_eval = self.xU(:, ind);
                self.xU(:, ind) = [];
            else
                while 1
                    constantSimplexSearch(self);
                    xm_grid = round(self.xm * self.mesh_size) / self.mesh_size;
                    [~, success] = points_neighbers_find(xm_grid, self);
                    sd_xc_grid = self.inter_par.interpolate_eval(xm_grid) - self.K * self.yRange * mindis(xm_grid, self.xE);
                    if success == 0
                        self.xU = [self.xU, xm_grid];
                        self.yU = [self.yU; sd_xc_grid];
                    else
                        if mindis(xm_grid, self.xE) < 1e-6
                            self.refine_trigger = 1;
                            break
                        else    
                            if (size(self.yU, 1) > 0 && sd_xc_grid < min(self.yU)) || size(self.yU, 1) == 0
                                self.xm_eval = xm_grid;
                                break
                            else
                                % One of the support points has the lower
                                % discrete search, add the newly found point to
                                % support points, and evaluated the minima of
                                % support points
                                self.xU = [self.xU, xm_grid];
                                [~, ind] = min(self.yU);
                                self.xm_eval = self.xU(:, ind);
                                self.xU(:, ind) = [];
                                break
                            end
                        end
                    end
                end
            end
        end
        
        function constantSimplexSearch(self)
            % Minimize the constant surrogate model
            self.xall = [self.xE, self.xU];
            self.tri = delaunayn(self.xall');
            Sc = zeros(1, size(self.tri, 1));
            Scl = zeros(1, size(self.tri, 1));
            for ii = 1 : size(self.tri, 1)
                [xc, R2] = circhyp(self.xall(:, self.tri(ii, :)), self.n);
                if R2 ~= inf
                    % initialization with body center of each simplex
                    center = self.xall(:, self.tri(ii,:)) * ones(self.n + 1, 1) / (self.n + 1); 
                    Sc(ii) = self.inter_par.interpolate_eval(center) - self.K * self.yRange * (R2-norm(center - xc)^2);
                    Scl(ii) = Sc(ii);
                    if ismember(self.ind_min, self.tri(ii, :)) ~= 1
                        Scl(ii) = inf;
                    end
                else
                    Sc(ii) = inf;
                    Scl(ii) = inf;
                end                    
            end
            [~, ind] = min(Sc);
            [x1, y1] = constantSearch(self, ind);

            [~, ind] = min(Scl);
            [x2, y2] = constantSearch(self, ind);
            if y2 < 2 * y1
                self.xm = x2;
            else
                self.xm = x1;
            end
        end
        
        function [x, y] = constantSearch(self, ind)
            % constantSearch: Handle function which minimizes the adaptive
            % surrogaet model with the given Delauany simplex, indexed as 
            % 'ind'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ind       :  The index of the Delauany simplex of interest
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Return
            % x         :  The minimizer of the surrogate model in this
            %              Delauany simplex
            % y         :  The surrogate model value at x
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [xc,R2] = circhyp(self.xall(:, self.tri(ind,:)), self.n);
            x0 = self.xall(:, self.tri(ind,:)) * ones(self.n + 1, 1) / (self.n + 1);
            fun = @(x) constantCost(self, R2, xc, x);
            options = optimoptions(@fmincon,'Algorithm','sqp','GradObj','On','DerivativeCheck','Off', 'Display', 'off');
            [x, y] = fmincon(fun, x0, [], [], [], [], x0*0, x0*0+1,[], options);
        end
        
        function [s, ds] = constantCost(self, R2, xc, x)
            % constantCost: Handle function which evaluates the cost of
            % constant K model at given x.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R2        :  The square of the radius of the circumsphere
            % xc        :  The circumcenter of the circumsphere
            % x         :  Given query point x of interest
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Return
            % s         :  The value of the Constant K surrogate model at x
            % ds        :  The derivative of the Constant K surrogate model
            %              at x
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            e = R2 - (x - xc)' * (x - xc);
            s = self.inter_par.interpolate_eval(x) - self.K * self.yRange * e;
            ds = self.inter_par.interpolate_grad(x) + 2 * self.K * self.yRange * (x - xc);
        end
        
        function adaptiveDelaunayMinimization(self)
            self.yU = zeros(size(self.xU, 2), 1);
            for ii = 1 : size(self.xU, 2)
                self.yU(ii) = (self.inter_par.interpolate_eval(self.xU(:, ii)) - self.y0) / mindis(self.xU(:, ii), self.xE);
            end
            if size(self.xU, 2) > 0 && min(self.yU) < 0
                [~, ind] = min(self.yU);
                self.xm_eval = self.xU(:, ind);
                self.xU(:, ind) = [];
            else
                while 1
                    adaptiveSimplexSearch(self);
                    xm_grid = round(self.xm * self.mesh_size) / self.mesh_size;
                    [~, success] = points_neighbers_find(xm_grid, self);
                    sd_xc_grid = (self.inter_par.interpolate_eval(xm_grid) - self.y0) / mindis(xm_grid, self.xE);
                    if success == 0
                        self.xU = [self.xU, xm_grid];
                        self.yU = [self.yU; sd_xc_grid];
                    else
                        if mindis(xm_grid, self.xE) < 1e-6
                            self.refine_trigger = 1;
                            break
                        else
                            if (size(self.yU, 1) > 0 && sd_xc_grid < min(self.yU)) || size(self.yU, 1) == 0
                                self.xm_eval = xm_grid;
                                break
                            else
                                % One of the support points has the lower
                                % discrete search, add the newly found point to
                                % support points, and evaluated the minima of
                                % support points
                                self.xU = [self.xU, xm_grid];
                                [~, ind] = min(self.yU);
                                self.xm_eval = self.xU(:, ind);
                                self.xU(:, ind) = [];
                                break
                            end
                        end
                    end
                end
            end
        end
        function adaptiveSimplexSearch(self)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the next point to evaluate through minimizing the 
            % adaptive surrogate model built upon Delaunay triangulation,
            % including evaluate the support points on the constraints.
            
            % adaptive surrogate model: s(x) = (p(x) - y0) / e(x)
            % Input :  Evaluated data points, xE and yE
            % Return:  Next point to evaluate, xm_eval
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            self.xall = [self.xE, self.xU];
            self.tri = delaunayn(self.xall');
            Sc = zeros(1, size(self.tri, 1));
            Scl = zeros(1, size(self.tri, 1));
            for ii = 1 : size(self.tri, 1)
                [xc, R2] = circhyp(self.xall(:, self.tri(ii, :)), self.n);
                if R2 ~= inf
                    % initialization with body center of each simplex
                    center = self.xall(:, self.tri(ii,:)) * ones(self.n + 1, 1) / (self.n + 1); 
                    Sc(ii) = (self.inter_par.interpolate_eval(center) - self.y0) / (R2-norm(center - xc)^2);
                    Scl(ii) = Sc(ii);
                    if ismember(self.ind_min, self.tri(ii, :)) ~= 1
                        Scl(ii) = inf;
                    end
                else
                    Sc(ii) = inf;
                    Scl(ii) = inf;
                end                    
            end
            [~, ind] = min(Sc);
            [x1, y1] = adaptiveSearch(self, ind);
            
            [~, ind] = min(Scl);
            [x2, y2] = adaptiveSearch(self, ind);
            if y2 < 2 * y1
                self.xm = x2;
                self.ym = y2;
            else
                self.xm = x1;
                self.ym = y1;
            end
        end
        function [x, y] = adaptiveSearch(self, ind)
            % adaptiveSearch: Handle function which minimizes the adaptive
            % surrogaet model with the given Delauany simplex, indexed as 
            % 'ind'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ind       :  The index of the Delauany simplex of interest
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Return
            % x         :  The minimizer of the surrogate model in this
            %              Delauany simplex
            % y         :  The surrogate model value at x
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [xc,R2] = circhyp(self.xall(:, self.tri(ind,:)), self.n);
            x0 = self.xall(:, self.tri(ind,:)) * ones(self.n + 1, 1) / (self.n + 1);
            fun = @(x) adaptiveCost(self, R2, xc, x);
            options = optimoptions(@fmincon,'Algorithm','sqp','GradObj','On','DerivativeCheck','Off', 'Display', 'off');
            [x, y] = fmincon(fun, x0, [], [], [], [], x0*0, x0*0+1,[], options);
        end
        function [s, ds] = adaptiveCost(self, R2, xc, x)
            % adaptiveCost: Handle function which evaluates the cost of
            % adaptive K model at given x.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R2        :  The square of the radius of the circumsphere
            % xc        :  The circumcenter of the circumsphere
            % x         :  Given query point x of interest
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Return
            % s         :  The value of the adaptive K surrogate model at x
            % ds        :  The derivative of the adaptive K surrogate model
            %              at x
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            e = R2 - (x - xc)' * (x - xc);
            p = self.inter_par.interpolate_eval(x);
            if p < self.y0
                s = -inf;
                ds = zeros(self.n, 1);
            else
                s = - e / (p - self.y0);
                dp = self.inter_par.interpolate_grad(x);
                de = -2 * (x - xc);
                ds = -de / (p - self.y0) + e * dp / (p - self.y0)^2;
            end                
        end
        
%         function plot_maker2D(self)
%             self.obj_lim = [min(self.surrogate_values), max(self.surrogate_values)];
%             N = size(self.DT_centers, 1);
%             X = zeros(N, 1);
%             Y = zeros(N, 1);
%             Z = zeros(N, 1);
%             for i = 1 : N
%                 point = self.DT_centers(i, :)';
%                 X(i) = point(1);
%                 Y(i) = point(2);
%                 Z(i) = self.surrogate_values(i);
%             end
%             figure; 
%             
%             % scatter plot the center of each Delaunay triangulation
%             scatter3(X, Y, Z, 30, 'filled')
%             
%             zlim(self.obj_lim);
%             hold on; grid on;
%             [~, ind_min] = min(Z);
%             
%             % scatter plot the minimizer of all centers
%             scatter3(X(ind_min), Y(ind_min), Z(ind_min), 'r', 'filled')
%             scatter3(X, Y, self.obj_lim(1) * ones(N, 1), 15, 'filled')
%             
%             % plot the vertical line of each center, connecting from the
%             % search model function values to the obj_lim lower bound
%             for i = 1 : N
%                 line([X(i); X(i)], [Y(i); Y(i)], [Z(i), self.obj_lim(1)], 'Color', 'black', 'Linestyle', '--');
%             end
%             
%             % plot the global minimum
%             scatter3(self.xmin(1), self.xmin(2), self.obj_lim(1), 'r*')
%             
%             % plot the quantized point of minimizer of all centers
%             scatter3(self.xc_eval(1), self.xc_eval(2), self.obj_lim(1), 'g', 'filled')
%             view(45, 45)
%             DelaunayTriangulationPlot2D(self);
%             if self.mesh_size > 16
%                 set(gca,'XTickLabel',[])
%                 set(gca,'YTickLabel',[])
%             end
%             set(gca,'ZTickLabel',[])
%             set(gca, 'xtick', [0:(1/self.mesh_size):1])
%             set(gca, 'ytick', [0:(1/self.mesh_size):1])
%             saveas(gca, ['figures/2D_DT_', num2str(self.iter), '.png'])
%         end
%         function DelaunayTriangulationPlot2D(self)
%             % plot the boundary lines of each Delaunay triangulation
%             N = size(self.DT_centers, 1);
%             comb = nchoosek(1:3, 2);
%             for i = 1 : N
%                 simplex = self.xE(:, self.tri(i, :));
%                 for j = 1 : size(comb, 1)
%                     line(simplex(1, comb(j, :)), simplex(2, comb(j, :)), [self.obj_lim(1), self.obj_lim(1)], 'Color', 'green');
%                 end
%             end
%         end
        function DelaunayUpdate(self)
            
            % update the evaluated point at the end of each iteration
            if self.refine_trigger ~= 1
                self.xE = [self.xE, self.xm_eval]; 
                physical_xm_eval = recover_physical_bounds(self.xm_eval, self.physical_lb, self.physical_ub);
                self.yE = [self.yE, self.func_eval(physical_xm_eval)];
            else
                self.mesh_size = self.mesh_size * 2;
                % if using constant K surrogate
                self.K = self.K * 2;
            end
        end
        function IterOptInfoOutput(self)
            fprintf('================== Iteration \t %i ================== \n', self.num_iter)
            fprintf('Mesh refine trigger \t %i \n', self.refine_trigger)
            fmt = ['Point to eval: [', repmat('%g, ', 1, numel(self.xm_eval)-1), '%g]\n'];
            fprintf(fmt, self.xm_eval)
            fprintf('Mesh size \t %i \n', self.mesh_size)
            fmt2 = ['Best point: [', repmat('%g, ', 1, numel(self.xE(:, find(self.yE==min(self.yE))))-1), '%g]\n'];
            fprintf(fmt2, self.xE(:, find(self.yE==min(self.yE))));
        end
    end
end





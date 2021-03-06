function [newadd, success] = points_neighbers_find(x, self)
    % Add the new point to the set
    % xE: evaluation points.
    % xU: unevaluated points.
    %
    %Author: Shahrouz Alimohammadi 2016, modified by Muhan 2021

    n = self.n;
    Ain = self.Ain;
    bin = self.bin;

    %keyboard
    % Find closest point to x
    [~, ~, x1]=mindis(x,[self.xE, self.xU]);


    % Caclulate the active constraints at x and x1
    %keyboard
    active_cons=1:2*n; b=bin-Ain*x; active_cons=active_cons(b<1e-3); 
    active_cons1=1:2*n; b=bin-Ain*x1; active_cons1=active_cons1(b<1e-3); 

    % Check the closest point is acceptable 
       if (isempty(active_cons) || min(ismember(active_cons,active_cons1))==1)
        % Acceptable case
               newadd=1;
               success=1;
               if mindis(x, self.xU)==0
                   newadd=0;
               end
        % New point is close to xU, and x1 is closer to xE. 
           %  if index>size(xE)
           %    newadd=0; x=x1;
           %    xU(:,index-size(xE))=[];
           %  end
       else
        % unacceptable case
        success = 0;
        newadd = 0; 
       end
end
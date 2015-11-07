function soln = hermiteSimpson(problem)
% soln = hermiteSimpson(problem)
%
% This function transcribes a trajectory optimization problem using the
% Hermite-Simpson (Seperated) method for enforcing the dynamics. It can be
% found in chapter four of Bett's book:
%
%   John T. Betts, 2001
%   Practical Methods for Optimal Control Using Nonlinear Programming
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'hermiteSimpson'
%   problem.options.hermiteSimpson = struct with method parameters:
%       .nGrid = number of grid points to use for transcription
%
% This transcription method is compatable with analytic gradients. To
% enable this option, set:
%   problem.nlpOpt.GradObj = 'on'
%   problem.nlpOpt.GradConstr = 'on'
%
% Then the user-provided functions must provide gradients. The modified
% function templates are as follows:
%
%         [dx, dxGrad] = dynamics(t,x,u)
%                 dx = [nState, nTime] = dx/dt = derivative of state wrt time
%                 dxGrad = [nState, 1+nx+nu, nTime]
%
%         [dObj, dObjGrad] = pathObj(t,x,u)
%                 dObj = [1, nTime] = integrand from the cost function
%                 dObjGrad = [1+nx+nu, nTime]
%
%         [c, ceq, cGrad, ceqGrad] = pathCst(t,x,u)
%                 c = [nCst, nTime] = column vector of inequality constraints  ( c <= 0 )
%                 ceq = [nCstEq, nTime] = column vector of equality constraints ( c == 0 )
%                 cGrad = [nCst, 1+nx+nu, nTime];
%                 ceqGrad = [nCstEq, 1+nx+nu, nTime];
%
%         [obj, objGrad] = bndObj(t0,x0,tF,xF)
%                 obj = scalar = objective function for boundry points
%                 objGrad = [1+nx+1+nx, 1]
%
%         [c, ceq, cGrad, ceqGrad] = bndCst(t0,x0,tF,xF)
%                 c = [nCst,1] = column vector of inequality constraints  ( c <= 0 )
%                 ceq = [nCstEq,1] = column vector of equality constraints ( c == 0 )
%                 cGrad = [nCst, 1+nx+1+nx];
%                 ceqGrad = [nCstEq, 1+nx+1+nx];
%

% Each segment needs an additional data point in the middle, thus:
nGrid = 2*problem.options.hermiteSimpson.nSegment-1;

% Print out some solver info if desired:
if problem.options.verbose > 0
    fprintf('  -> Transcription via Hermite-Simpson method, nSegment = %d\n',...
        problem.options.hermiteSimpson.nSegment);
end

%%%% Method-specific details to pass along to solver:

% Hermite-Simpson calculation of defects:
problem.func.defectCst = @computeDefects;

%%%% The key line - solve the problem by direct collocation:
soln = directCollocation(problem);

% Use piecewise cubic interpolation for each trajectory segment
tSoln = soln.grid.time';
xSoln = soln.grid.state';
uSoln = soln.grid.control';
soln.interp.state = @(t)( interp1(tSoln,xSoln,t,'pchip')' );
soln.interp.control = @(t)( interp1(tSoln,uSoln,t,'pchip')' );

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                          SUB FUNCTIONS                            %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%



function [defects, defectsGrad] = computeDefects(dt,x,f,dtGrad,xGrad,fGrad)
%
% This function computes the defects that are used to enforce the
% continuous dynamics of the system along the trajectory.
%
% INPUTS:
%   dt = time step (scalar)
%   x = [nState, nTime] = state at each grid-point along the trajectory
%   f = [nState, nTime] = dynamics of the state along the trajectory
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   dtGrad = [2,1] = gradient of time step with respect to [t0; tF]
%   xGrad = [nState,nTime,nDecVar] = gradient of trajectory wrt dec vars
%   fGrad = [nState,nTime,nDecVar] = gradient of dynamics wrt dec vars
%
% OUTPUTS:
%   defects = [nState, nTime-1] = error in dynamics along the trajectory
%   defectsGrad = [nState, nTime-1, nDecVars] = gradient of defects
%

nTime = size(x,2);

iLow = 1:2:(nTime-1);
iMid = iLow + 1;
iUpp = iMid + 1;

xLow = x(:,iLow);
xMid = x(:,iMid);
xUpp = x(:,iUpp);

fLow = f(:,iLow);
fMid = f(:,iMid);
fUpp = f(:,iUpp);

% Mid-point constraint (Hermite)
defectMidpoint = xMid - (xUpp+xLow)/2 - dt*(fLow-fUpp)/4;

% Interval constraint (Simpson)
defectInterval = xUpp - xLow - dt*(fUpp + 4*fMid + fLow)/3;

% Pack up all defects:
defects = [defectMidpoint, defectInterval];


%%%% Gradient Calculations:
if nargout == 2
    
    xLowGrad = xGrad(:,iLow,:);
    xMidGrad = xGrad(:,iMid,:);
    xUppGrad = xGrad(:,iUpp,:);
    
    fLowGrad = fGrad(:,iLow,:);
    fMidGrad = fGrad(:,iMid,:);
    fUppGrad = fGrad(:,iUpp,:);
    
    % Mid-point constraint (Hermite)
    dtGradTerm = zeros(size(xMidGrad));
    dtGradTerm(:,:,1) = -dtGrad(1)*(fLow-fUpp)/4;
    dtGradTerm(:,:,2) = -dtGrad(2)*(fLow-fUpp)/4;
    defectMidpointGrad = xMidGrad - (xUppGrad+xLowGrad)/2 + dtGradTerm + ...
        - dt*(fLowGrad-fUppGrad)/4;
    
    % Interval constraint (Simpson)
    dtGradTerm = zeros(size(xUppGrad));
    dtGradTerm(:,:,1) = -dtGrad(1)*(fUpp + 4*fMid + fLow)/3;
    dtGradTerm(:,:,2) = -dtGrad(2)*(fUpp + 4*fMid + fLow)/3;
    defectIntervalGrad = xUppGrad - xLowGrad + dtGradTerm + ...
        - dt*(fUppGrad + 4*fMidGrad + fLowGrad)/3;
    
    %Pack up the gradients of the defects:
    defectsGrad = cat(2,defectMidpointGrad,defectIntervalGrad);
    
end

end


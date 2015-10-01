function soln = multipleShooting(problem)
% soln = multipleShooting(problem)
%
% TrajOpt utility function. Implements the core of the multiple shooting
% algorithms, including 4th-order runge-kutta and also the mid-point method
%

%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

% Figure out grid size:
nSegment = Opt.multipleShooting.nSegment;
nSubStep = Opt.multipleShooting.nSubStep;
nGridControl = 2*nSegment*nSubStep + 1;
nGridState = nSegment + 1;

% Print out some solver info if desired:
if Opt.verbose > 0
    fprintf('        nSegments = %d \n', nSegment);
    fprintf('        nSubSteps = %d \n', nSubStep);
end

% Interpolate the guess at the transcription grid points for initial guess:
guess.tSpan = G.time([1,end]);
guess.tState = linspace(guess.tSpan(1), guess.tSpan(2), nGridState);
guess.tControl = linspace(guess.tSpan(1), guess.tSpan(2), nGridControl);
guess.state = interp1(G.time', G.state', guess.tState')';
guess.control = interp1(G.time', G.control', guess.tControl')';
[zGuess, pack] = packDecVar(guess.tSpan, guess.state, guess.control);

% Unpack all bounds:
tLow = [B.initialTime.low, B.finalTime.low];
xLow = [B.initialState.low, B.state.low*ones(1,nGridState-2), B.finalState.low];
uLow = B.control.low*ones(1,nGridControl);
zLow = packDecVar(tLow,xLow,uLow);

tUpp = [B.initialTime.upp, B.finalTime.upp];
xUpp = [B.initialState.upp, B.state.upp*ones(1,nGridState-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nGridControl);
zUpp = packDecVar(tUpp,xUpp,uUpp);

%%%% Set up problem for fmincon:

dynFun = @(t,x,u)( combinedDynamics(t,x,u,F.dynamics,F.pathObj) );
intFun = problem.options.multipleShooting.integralFun;

integralFun = @(dt, t0, x0, uLow, uMid, uUpp)intFun( dt, t0, x0, uLow, uMid, uUpp, dynFun);

P.objective = @(z)( ...
    myObjective(z, pack, integralFun, F.bndObj) );

P.nonlcon = @(z)( myConstraint(z, pack, integralFun, F.pathCst, F.bndCst) );

P.x0 = zGuess;
P.lb = zLow;
P.ub = zUpp;
P.Aineq = []; P.bineq = [];
P.Aeq = []; P.beq = [];
P.options = Opt.nlpOpt;
P.solver = 'fmincon';

%%%% Call fmincon to solve the non-linear program (NLP)
tic;
[zSoln, objVal,exitFlag,output] = fmincon(P);
[tSpan,~,uSoln] = unPackDecVar(zSoln,pack);
nlpTime = toc;

%%%% Store the results:
[tGrid,xGrid,uGrid] = simulateSystem(zSoln, pack, integralFun);
soln.grid.time = tGrid;
soln.grid.state = xGrid;
soln.grid.control = uGrid;

% Quadratic interpolation over each sub-step for the control:
tSoln = linspace(tSpan(1),tSpan(2),nGridControl);
soln.interp.control = @(t)( interp1(tSoln', uSoln', t','pchip')' );

% Cubic spline representation of the state over each substep:
dxGrid = F.dynamics(tGrid,xGrid,uGrid);
xSpline = pwch(tGrid, xGrid, dxGrid);
soln.interp.state = @(t)( ppval(xSpline,t) );

% General information about the optimization run
soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                   SUB FUNCTIONS                                   %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [decVars,pack] = packDecVar(tSpan,state,control)
%
% This function collapses the time (t), state (x)
% and control (u) matricies into a single vector
%
% INPUTS:
%   tSpan = [1, 2] = time bounds
%   state = [nState, nGridState] = state vector at each grid point
%   control = [nControl, nGridControl] = control vector at each grid point
%
% OUTPUTS:
%   decVars = column vector of 2 + nState*nGridState + nControl*nGridControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nState
%       .nGridState
%       .nControl
%       .nGridControl
%
% NOTES:
% nGridControl = 2*nSegment*nSubStep + 1;
% nGridState = nSegment + 1;
%

[nState, nGridState] = size(state);
[nControl, nGridControl] = size(control);

xCol = reshape(state, nState*nGridState, 1);
uCol = reshape(control, nControl*nGridControl, 1);

decVars = [tSpan(:);xCol;uCol];

pack.nState = nState;
pack.nGridState = nGridState;
pack.nControl = nControl;
pack.nGridControl = nGridControl;
pack.nSegment = nGridState - 1;
pack.nSubStep = (nGridControl-1)/(2*pack.nSegment);

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [tSpan, state, control] = unPackDecVar(decVars,pack)
%
% This function unpacks the decision variables for
% trajectory optimization into the time (t),
% state (x), and control (u) matricies
%
% INPUTS:
%   decVars = column vector of 2 + nState*nGridState + nControl*nGridControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nState
%       .nGridState
%       .nControl
%       .nGridControl
%
% OUTPUTS:
%   tSpan = [1, 2] = time bounds
%   state = [nState, nGridState] = state vector at each grid point
%   control = [nControl, nGridControl] = control vector at each grid point
%

nx = pack.nState*pack.nGridState;
nu = pack.nControl*pack.nGridControl;

tSpan = [decVars(1),decVars(2)];

state = reshape(decVars((2+1):(2+nx)), pack.nState, pack.nGridState);
control = reshape(decVars((2+nx+1):(2+nx+nu)), pack.nControl, pack.nGridControl);

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function cost = myObjective(decVars, pack,integralFun, bndObj)
%
% This function unpacks the decision variables, sends them to the
% user-defined objective functions, and then returns the final cost
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%   bndObj = user-defined boundary objective function
%
% OUTPUTS:
%   cost = scalar cost for this set of decision variables
% 
%

% All of the real work happens inside this function:
[t,x,~,~,pathCost] = simulateSystem(decVars, pack, integralFun);

% Compute the cost at the boundaries of the trajectory
if isempty(bndObj)
    bndCost = 0;
else
    t0 = t(1);
    tF = t(end);
    x0 = x(:,1);
    xF = x(:,end);
    bndCost = bndObj(t0,x0,tF,xF);
end

cost = bndCost + pathCost;

end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [c, ceq] = myConstraint(decVars, pack, integralFun, pathCst, bndCst)
%
% This function unpacks the decision variables, computes the defects along
% the trajectory, and then evaluates the user-defined constraint functions.
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%   pathCst = user-defined path-constraint function
%   bndCst = user-defined boundary constraint function
%
% OUTPUTS:
%   c = non-linear inequality constraint
%   ceq = non-linear equatlity cosntraint
%
% NOTE:
%   - path constraints are  satisfied at the start and end of each sub-step
%


[t,x,u,defects] = simulateSystem(decVars, pack, integralFun);

%%%% Call user-defined constraints and pack up:
[c, ceq] = collectConstraints(t,x,u,...
    defects,...
    pathCst, bndCst);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [t,x,u,defects,pathCost] = simulateSystem(decVars, pack, integralFun)
%
% This function does the real work of the transcription method. It
% simulates the system forward in time across each segment of the
% trajectory, computes the integral of the cost function, and then matches
% up the defects between the end of each segment and the start of the next.
%
% INPUTS:
%   decVars = column vector of decision variables
%   pack = details about how to convert decision variables into t,x, and u
%   dynamics = user-defined dynamics function handle
%   pathObj = user-defined path-objective function
%
% OUTPUTS:
%   t = [1 x nGrid] = time vector for the edges of the sub-step grid
%   x = [nState x nGrid] = state vector
%   u = [nControl x nGrid] = control vector
%   defects = [nState x nSegment] = defect matrix
%   pathCost = scalar cost for the path integral
%   
% NOTES:
%   - nGrid = nSegment*nSubStep+1
%   - This function is usually called twice for each combination of
%   decision variables: once by the objective function and once by the
%   constraint function. To keep the code fast I cache the old values and
%   only recompute when the inputs change.
%   


%%%% CODE OPTIMIZATION %%%%
%
% Prevents the same exact code from being called twice by caching the
% solution and reusing it when appropriate.
%
global MULT_SHOOT_t MULT_SHOOT_x MULT_SHOOT_u
global MULT_SHOOT_defects MULT_SHOOT_pathCost
global MULT_SHOOT_decVars
%
usePreviousValues = false;
if ~isempty(MULT_SHOOT_decVars)
    if length(MULT_SHOOT_decVars) == length(decVars)
        if ~any(MULT_SHOOT_decVars ~= decVars)
            usePreviousValues = true;
        end
    end
end
%
if usePreviousValues
    t = MULT_SHOOT_t;
    x = MULT_SHOOT_x;
    u = MULT_SHOOT_u;
    defects = MULT_SHOOT_defects;
    pathCost = MULT_SHOOT_pathCost;
else
%
%
%%%% END CODE OPTIMIZATION %%%%


    [tSpan, state, control] = unPackDecVar(decVars,pack);
    
    nState = pack.nState;
    nSegment = pack.nSegment;
    nSubStep = pack.nSubStep;
    
    % NOTES:
    %   The following bit of code is a bit confusing, mostly due to the
    %   need for vectorization to make things run at a reasonable speed in
    %   Matlab. Part of the confusion comes because the decision variables
    %   include the state at the beginning of each segment, but the control
    %   at the beginning and middle of each substep - thus there are more
    %   control grid-points than state grid points. The calculations are
    %   vectorized over segments, but not sub-steps, since the result of
    %   one sub-step is required for the next.
    
    % time, state, and control at the ends of each substep
    nTime = 1+nSegment*nSubStep;
    t = linspace(tSpan(1), tSpan(2), nTime);
    x = zeros(nState, nTime);
    u = control(:,1:2:end); % Control a the endpoints of each segment
    uMid = control(:,2:2:end);  %Control at the mid-points of each segment
    c = zeros(1, nTime-1);  %Integral cost for each segment
    dt = (t(end)-t(1))/(nTime-1);
       
    idx = 1:nSubStep:(nTime-1);   %Indicies for the start of each segment
    x(:,[idx,end]) = state;   %Fill in the states that we already know
    for iSubStep = 1:nSubStep
        % March forward Runge-Kutta step    
        
        t0 = t(idx);
        x0 = x(:,idx);
        
        [xNext, c(idx)] = integralFun(dt, t0, x0, u(:,idx), uMid(:,idx), u(:,idx+1));
                
        if iSubStep == nSubStep %We've reached the end of the interval
            % Compute the defect vector:
            defects = xNext - x(:,idx+1);
        else
            % Store the state for next step in time
            idx = idx+1;   %  <-- This is important!!
            x(:,idx) = xNext;
        end
        
    end
    
    pathCost = sum(c);  %Sum up the integral cost over each segment
    
    %%%% Cache results to use on the next call to this function.
    MULT_SHOOT_t = t;
    MULT_SHOOT_x = x;
    MULT_SHOOT_u = u;
    MULT_SHOOT_defects = defects;
    MULT_SHOOT_pathCost = pathCost;
    
end

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function dz = combinedDynamics(t,x,u,dynamics,pathObj)
% dz = combinedDynamics(t,x,u,dynamics,pathObj)
%
% This function packages the dynamics and the cost function together so
% that they can be integrated at the same time.
%
% INPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%   dynamics(t,x,u) = dynamics function handle
%               dx = [nState, nTime] = dx/dt = derivative of state wrt time
%   pathObj(t,x,u) = integral cost function handle
%                 dObj = [1, nTime] = integrand from the cost function
%
% OUTPUTS:
%   dz = [dx; dObj] = combined dynamics of state and cost
%

dx = dynamics(t,x,u);
if isempty(pathObj)
    dc = zeros(size(t));
else
    dc = pathObj(t,x,u);
end

dz = [dx;dc];  %Combine and return

end



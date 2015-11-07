function soln = radauCollocation(problem)
% soln = radauCollocation(problem)
%
% This function transcribes a trajectory optimization problem using
% Radau Orthogonal Collocation, as described by Rao et al:
%
% "A unified framework for the numerical solution of optimal control
% problems using pseudospectral methods". Divya Garg, Michael Patterson,
% WIlliam W. Hagar, Anil V. Rad, David A. Benson, Geoffrey T. Huntington.
% Automatica 2010, Elsevier
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'radauCollocation'
%   problem.options.radauCollocation = struct with method parameters:
%       .nColPts = number of collocation points in the transcription
%

%To make code more readable
G = problem.guess;
B = problem.bounds;
F = problem.func;
Opt = problem.options;

%Internal Parameters:
P.n = Opt.radauCollocation.nColPts;
P.n1 = P.n+1;
P.x = [radaupts(P.n); 1];   % Barycentric points
[P.D, P.v] = collocD(P.x);   %Full differentation matrix and barycentric weights
P.Dc = P.D(1:(end-1),:);  %Differentiation matrix for collocation points

nGrid = P.n1;

if isempty(F.pathObj)
    error('Empty path integral not yet supported in radauCollocation');
else
    P.dynFun = @(t,y,u)( dynFun(t,y,u, F.pathObj, F.dynamics) );
    
    objGuess = zeros(1,nGrid);
    guess.tSpan = G.time([1,end]);
    guess.time = xScale(P.x,guess.tSpan)';
    guess.state = [objGuess;     % prepend objective function to state
        interp1(G.time', G.state', guess.time')'];
    guess.control = interp1(G.time', G.control', guess.time')';
    
    
    [zGuess, pack] = packDecVar(guess.time, guess.state, guess.control);
    
    % Unpack all bounds:
    dummyMatrix = zeros(1,nGrid-2);  %This just needs to be the right size
    
    objLow = [0, -inf(1,nGrid-1)];
    tLow = [B.initialTime.low, dummyMatrix, B.finalTime.low];
    xLow = [objLow;
        B.initialState.low, B.state.low*ones(1,nGrid-2), B.finalState.low];
    uLow = B.control.low*ones(1,nGrid);
    zLow = packDecVar(tLow,xLow,uLow);
    
    objUpp = [0, inf(1,nGrid-1)];
    tUpp = [B.initialTime.upp, dummyMatrix, B.finalTime.upp];
    xUpp = [objUpp;
        B.initialState.upp, B.state.upp*ones(1,nGrid-2), B.finalState.upp];
    uUpp = B.control.upp*ones(1,nGrid);
    zUpp = packDecVar(tUpp,xUpp,uUpp);
    
end


%%%% Set up problem for fmincon:

Problem.objective = @(z)( ...
    myObjective(z, pack, F.bndObj, P) );

Problem.nonlcon = @(z)( ...
    myConstraint(z, pack, F.pathCst, F.bndCst, P) );

Problem.x0 = zGuess;
Problem.lb = zLow;
Problem.ub = zUpp;
Problem.Aineq = []; P.bineq = [];
Problem.Aeq = []; P.beq = [];
Problem.options = Opt.nlpOpt;
Problem.solver = 'fmincon';


%%%% Call fmincon to solve the non-linear program (NLP)
tic;
[zSoln, objVal,exitFlag,output] = fmincon(Problem);
[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack,P);
nlpTime = toc;

%%%% Store the results:
soln.grid.time = tSoln;
soln.grid.state = xSoln(2:end,:);
soln.grid.control = uSoln;

%%%% Interpolate the results:
xxSoln = xScale(P.x,tSoln([1,end]));
soln.interp.state = @(t)( barycentricInterpolate(t', xSoln(2:end,:)',xxSoln,P.v)' );
soln.interp.control = @(t)( barycentricInterpolate(t', uSoln',xxSoln,P.v)' );

%%%% Solution information:
soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;

soln.problem = problem;  % Return the fully detailed problem struct
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                          Sub-Functions                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function dy = dynFun(t,y,u,pathObj, dynamics)
%
% This function is a wrapper for the user-defined dynamics and path
% integral functions

x = y(2:end,:); %First row is the cost integral

% Call user-defined functions:
dObj = pathObj(t,x,u);
dx = dynamics(t,x,u);

% Pack up:
dy = [dObj; dx];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = xScale(x,d)
%Rescale the points in x to an arbitrary interval, assuming that the input
%set of points are on the interval [-1,1];
x = 0.5*(d(2)+d(1) + (d(2)-d(1))*x);
end


function D = DScale(D,d)
%Rescale the differentiation matrix D to an arbitrary interval, assuming
% that the input set of points are on the interval [-1,1];
D = 2*D/(d(2)-d(1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z,pack] = packDecVar(t,y,u)
%
% This function collapses the time (t), state (y)
% and control (u) matricies into a single vector
%
% INPUTS:
%   t = [1, nTime] = time vector (grid points)
%   y = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%
% OUTPUTS:
%   z = column vector of 2 + nTime*(nState+nControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nTime
%       .nState
%       .nControl
%

nTime = length(t);
nState = size(y,1);
nControl = size(u,1);

tSpan = [t(1); t(end)];
yCol = reshape(y, nState*nTime, 1);
uCol = reshape(u, nControl*nTime, 1);

z = [tSpan;yCol;uCol];

pack.nTime = nTime;
pack.nState = nState;
pack.nControl = nControl;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [t,y,u] = unPackDecVar(z,pack,P)
%
% This function unpacks the decision variables for
% trajectory optimization into the time (t),
% state (x), and control (u) matricies
%
% INPUTS:
%   z = column vector of 2 + nTime*(nState+nControl) decision variables
%   pack = details about how to convert z back into t,x, and u
%       .nTime
%       .nState
%       .nControl
%
% OUTPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%

nTime = pack.nTime;
nState = pack.nState;
nControl = pack.nControl;
ny = nState*nTime;
nu = nControl*nTime;

t = xScale(P.x,z([1,2]));  %Rescale time
y = reshape(z((2+1):(2+ny)),nState,nTime);
u = reshape(z((2+ny+1):(2+ny+nu)),nControl,nTime);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cost = myObjective(z,pack,bndObj,P)

[t,y] = unPackDecVar(z,pack,P);

% The first row of the state (y) is the integral cost function
integralCost = y(1,end);


% Compute the cost at the boundaries of the trajectory
if isempty(bndObj)
    bndCost = 0;
else
    t0 = t(1);
    tF = t(end);
    x0 = y(2:end,1);
    xF = y(2:end,end);
    bndCost = bndObj(t0,x0,tF,xF);
end

cost = bndCost + integralCost;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [c, ceq] = myConstraint(z, pack, pathCst, bndCst, P)

[t,y,u] = unPackDecVar(z,pack,P);
tSpan = t([1,end]);
Dc = DScale(P.Dc,tSpan);

% Compute the derivative by differentiating the state:
dyFun = (Dc*(y'))';

% Compute the derivative by calling dynamics function:
idxCol = 1:(length(t)-1);  %Collocation points (ommit last point);
dyDyn = P.dynFun(t(idxCol),y(:,idxCol),u(:,idxCol));  

% Defects by matching the derivative at each collocation point:
defects = dyFun - dyDyn;

%%%% Call user-defined constraints and pack up:
x = y(2:end,:);  %First row is the cost function - ommit
[c, ceq] = collectConstraints(t,x,u,defects, pathCst, bndCst);

end


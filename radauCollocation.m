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
    P.dynFun = @(t,y,u)( dynFun(t,y,u,pathObj, dynamics) );
    
    objGuess = zeros(1,nGrid);
    guess.tSpan = G.time([1,end]);
    guess.time = xScale(P.x,guess.tSpan);
    guess.state = [objGuess;     % prepend objective function to state 
        interp1(G.time', G.state', guess.time')'];
    guess.control = interp1(G.time', G.control', guess.time')';    

% Unpack all bounds:
dummyMatrix = zeros(1,nColPts-2);  %This just needs to be the right size

objLow = -inf(1,nGrid);
tLow = [B.initialTime.low, dummyMatrix, B.finalTime.low];
xLow = [objLow;
    B.initialState.low, B.state.low*ones(1,nColPts-2), B.finalState.low];
uLow = B.control.low*ones(1,nColPts);
zLow = packDecVar(tLow,xLow,uLow);

objUpp = inf(1,nGrid);
tUpp = [B.initialTime.upp, dummyMatrix, B.finalTime.upp];
xUpp = [objUpp;
    B.initialState.upp, B.state.upp*ones(1,nColPts-2), B.finalState.upp];
uUpp = B.control.upp*ones(1,nColPts);
zUpp = packDecVar(tUpp,xUpp,uUpp);

end


%%%% Set up problem for fmincon:

P.objective = @(z)( ...
    myObjective(z, pack, F.bndObj, P) );

P.nonlcon = @(z)( ...
    myConstraint(z, pack, F.pathCst, F.bndCst, P) );

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
[tSoln,xSoln,uSoln] = unPackDecVar(zSoln,pack,P);
nlpTime = toc;

%%%% Store the results:
soln.grid.time = tSoln;
soln.grid.state = xSoln(2:end,:);
soln.grid.control = uSoln;

%%%% Interpolate the results:

%%%% TODO !!

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
x = 0.5*(d(2)+d(1) + (d(2)-d(1)))*x;
end


function D = DScale(D,x)
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

function obj = myObjective(z,pack,bndObj,P)


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, ceq] = myConstraint(z, pack, pathCst, bndCst, P)


end


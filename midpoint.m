function soln = midpoint(problem)
% soln = midpoint(problem)
%
% This function transcribes a trajectory optimization problem using the
% multiple shooting, using the "mid-point method" for integration over each
% individual time step.
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'midpoint'
%   problem.options.midpoint = struct with method parameters:
%       .nSegment = number of trajectory segments
%       .nSubStep = number of sub-steps to use in each segment
%

% Print out some solver info if desired:
if problem.options.verbose > 0
    fprintf('  -> Transcription via mid-point method \n');  
end

% Pack up things to send to the general-purpose multiple shooting solver:
problem.options.multipleShooting = problem.options.midpoint;
problem.options.multipleShooting.integralFun = @midpointCore;

soln = multipleShooting(problem);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                   SUB FUNCTIONS                                   %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [xNext, cost] = midpointCore(dt, t0, x0, uLow, uMid, ~, dynFun)
        
        nState = size(x0,1);

        k0 = dynFun(t0,        x0,                         uLow);
        k1 = dynFun(t0+0.5*dt, x0 + 0.5*dt*k0(1:nState,:), uMid);
        z = dt*k1;  %Change over the sub-step
          
        xNext = x0 + z(1:nState,:);  %Next state
        cost = z(end,:);  %Integral of the cost function over this step

end

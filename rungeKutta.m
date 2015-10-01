function soln = rungeKutta(problem)
% soln = rungeKutta(problem)
%
% This function transcribes a trajectory optimization problem using the
% multiple shooting, with 4th-order Runge Kutta integration
%
% See Bett's book for details on the method
%
% For details on the input and output, see the help file for trajOpt.m
%
% Method specific parameters:
%
%   problem.options.method = 'rungeKutta'
%   problem.options.rungeKutta = struct with method parameters:
%       .nSegment = number of trajectory segments
%       .nSubStep = number of sub-steps to use in each segment
%

% Print out some solver info if desired:
if problem.options.verbose > 0
    fprintf('  -> Transcription via 4th-order Runge-Kutta method \n');  
end

% Pack up things to send to the general-purpose multiple shooting solver:
problem.options.multipleShooting = problem.options.rungeKutta;
problem.options.multipleShooting.integralFun = @rungeKuttaCore;

soln = multipleShooting(problem);

end


%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%
%%%%                   SUB FUNCTIONS                                   %%%%
%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%


function [xNext, cost] = rungeKuttaCore(dt, t0, x0, uLow, uMid, uUpp, dynFun)
        
        nState = size(x0,1);

        k0 = dynFun(t0,        x0,                         uLow);
        k1 = dynFun(t0+0.5*dt, x0 + 0.5*dt*k0(1:nState,:), uMid);
        k2 = dynFun(t0+0.5*dt, x0 + 0.5*dt*k1(1:nState,:), uMid);
        k3 = dynFun(t0+dt,     x0 +     dt*k2(1:nState,:), uUpp);
        z = (dt/6)*(k0 + 2*k1 + 2*k2 + k3);  %Change over the sub-step
          
        xNext = x0 + z(1:nState,:);  %Next state
        cost = z(end,:);  %Integral of the cost function over this step

end

function y = barycentricInterpolate(x,yk,xk,vk)
% y = barycentricInterpolate(x,yk,xk,vk)
%
% Interpolates an orthogonal polynomial using barycentric interpolation
%
% INPUTS:
%   x = [nTime, 1] = vector of points to evaluate polynomial at
%   yk = [nGrid, nOut] = value of the function to be interpolated at each
%       grid point
%   xk = [nGrid, 1] = roots of orthogonal polynomial
%   vk = [nGrid, 1] = barycentric interpolation weights
%
% OUTPUTS:
%   y = [nTime, nOut] = value of the function at the desired points
%
% NOTES:
%   xk and yk should be produced by chebfun (help chebpts)
%

nOut = size(yk,2);
nTime = length(x);
y = zeros(nTime, nOut);

for i=1:nOut
    y(:,i) = bary(x,yk(:,i),xk,vk);
end

end




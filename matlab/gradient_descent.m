function [xi, yi] = gradient_descent(x, y, Z, xi, yi)
% Simple gradient descent routine for finding local minima.
%
% INPUT x : uniformly spaced x-axis (Nx elements)
%
% INPUT y : uniformly spaced y-axis (Ny elements)
%
% INPUT Z : a 2-D matrix of z values (a Ny-by-Nx matrix created using
% meshgrid).
%
% INPUT xi, yi : coordinate to interpolate Z at.
%
% OUTPUT xi, yi : position of local minima
%
% Author: John Gebbie
% Creation Date: 2013 Aug 7

assert(mean(abs(diff(x,2))) <= eps, 'x must be equally spaced');
assert(mean(abs(diff(y,2))) <= eps, 'y must be equally spaced');
assert(length(x) == size(Z,2));
assert(length(y) == size(Z,1));

[X,Y] = ndgrid(x,y);

dx = diff(x(1:2));
dy = diff(y(1:2));
[dZdy,dZdx] = gradient(Z.',dy,dx);

interp_dZdx = griddedInterpolant(X,Y,dZdx,'spline');
interp_dZdy = griddedInterpolant(X,Y,dZdy,'spline');

p = [xi yi];

N = 2*ceil(norm(size(Z)));

for n = 1:N
    
    dZdx1 = interp_dZdx(p);
    dZdy1 = interp_dZdy(p);
    
    u = [dZdx1*dx,dZdy1*dy];
    u = -u ./ norm(u);       % unit vector in direction of greatest descent
    
    alpha = (N-n+1)/N;
    
    p = p + [dx dy] .* u .* alpha;
end

xi = p(1);
yi = p(2);

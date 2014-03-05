function [yi] = smooth_triangle(x,y,xi,xhalfbase,no_norm)
% Simple triangle smoothing, interpolation, and extrapolation routine from
% irregularly sampled input data points to any output data points.
% 1-dimensional.
%
% INPUT x : x-coordinate of input data
%
% INPUT y : y-coordinate (value) of input data
%
% INPUT xi : x-coordinate to interpolate y-value at
%
% INPUT xhalfbase : half width of base of full triangle
%
% INPUT no_norm : do not normalize within a bin, sum instead
%
% OUTPUT yi : interpolated y-value. same dimensions as xi.
%
% Author: John Gebbie
% Modify Date: 2014 Mar 04

narginchk(4,5);
if nargin < 5
    no_norm = false;
end

[Xi, X] = ndgrid(xi(:), x(:));
D = abs(Xi - X);
V = -D./xhalfbase + 1;
V(V<0) = 0;
W = V;
if ~no_norm
    W = W ./ repmat(sum(W,2), 1, size(W,2));
end
W(~isfinite(W)) = 0;
yi = nan(size(xi));
yi(:) = W*y(:);

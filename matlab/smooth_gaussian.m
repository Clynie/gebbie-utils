function [yi] = smooth_gaussian(x,y,xi,xsigma)
% Simple Gaussian smoothing, interpolation, and extrapolation routine from
% irregularly smapled input data points to any output data points.
% 1-dimensional.
%
% INPUT x : x-coordinate of input data
%
% INPUT y : y-coordinate (value) of input data
%
% INPUT xi : x-coordinate to interpolate y-value at
%
% INPUT xsigma : sigma (width) of Gaussian
%
% OUTPUT yi : interpolated y-value. same dimensions as xi.
%
% Author: John Gebbie
% Modify Date: 2013 Oct 08

%normpdf = @(x,mu,sigma) exp(-(x-mu).^2./(2.*sigma.^2));
lognormpdf = @(x,mu,sigma) -(x-mu).^2./(2.*sigma.^2);

[Xi, X] = ndgrid(xi(:), x(:));
L = lognormpdf(Xi,X,xsigma);
L = L - repmat(max(L,[],2), 1, size(L,2));
V = exp(L);
W = V ./ repmat(sum(V,2), 1, size(V,2));
yi = nan(size(xi));
yi(:) = W*y(:);

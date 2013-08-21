function [fn_of_a] = deriv(f,n,a,delta)
% Nth derivative of a function, computed numerically
%
% INPUT f : function taking a single argument
%
% INPUT n : order of derivative
%
% INPUT a : range point at which to compute derivative
%
% INPUT delta : the separation from the range point to use to compute
% slopes
%
% OUTPUT fn_of_a : derivative at point
%
% Author: John Gebbie
% Modify Date: 2013 Aug 20

M = 5;
x = linspace(a-delta, a+delta, M);
y = nan(1,M);
for m = 1:M
    y(m) = f(x(m));
end
dx = diff(x(1:2));
for n1 = 1:n
    y = gradient(y, dx);
end
fn_of_a = y((M+1)/2);
end

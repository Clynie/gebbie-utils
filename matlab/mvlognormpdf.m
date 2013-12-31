function [v] = mvlognormpdf(x,mu,sigma)
% A multivariate normal PDF
%
% INPUT x - points to evaluate, a N-by-D matrix of N points each of D
% dimension
%
% INPUT mu - mean of multivariate normal, a 1-by-D vector
%
% INPUT sigma - covariance matrix, a D-by-D matrix
%
% OUTPUT v - output values, a N-by-1 vector.
%
% Author: John Gebbie
% Institution: Metron
% Creation Date: 2013, Dec 30

N = size(x,1);
D = size(x,2);

tmp = log( ( (2*pi)^D * det(sigma) )^(-1/2) );

v = nan(N,1);
for n = 1:N
    t = x(n,:) - mu;
    v(n) = tmp + -1/2 * t / sigma * t';
end

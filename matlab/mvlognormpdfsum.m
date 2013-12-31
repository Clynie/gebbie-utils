function [v] = mvlognormpdfsum(x,mu,sigma)
% Log sum of multivariate normal PDFs relative to a fixed mean.
%
% INPUT x - points to evaluate pdfs at, a N-by-D matrix of N points each of
% D dimension
%
% INPUT mu - mean of multivariate normal, a 1-by-D vector
%
% INPUT sigma - covariance matrix, a D-by-D matrix
%
% OUTPUT v - output value, a scalar
%
% Author: John Gebbie
% Institution: Metron
% Creation Date: 2013, Dec 30

t0 = mvlognormpdf(x,mu,sigma);
t1 = max(t0);
t2 = t0 - t1;
t3 = exp(t2);
t4 = sum(t3);
t5 = log(t4);
v = t5 + t1;

function [mu_x, Sigma_x, SP, SP_x] = unscented_transform(mu, Sigma, f)
% Unscented transform of an N-dimensional Guassian through an arbitrary
% function having a P-dimensional range.
%
% INPUT mu : mean (an N-by-1 vector)
%
% INPUT Sigma : covariance matrix (an N-by-N matrix)
%
% INPUT f : function taking M inputs points as an N-by-M matrix and
% outputing a P-by-M matrix of transformed points
%
% OUTPUT mu_x : transformed mean (a P-by-1 vector)
%
% OUTPUT Sigma_x : transformed covariance matrix (a P-by-P matrix)
%
% SP : sigma points (input domain)
%
% SP_x : sigma points (output domain)
%
% Author: John Gebbie
% Institution: Portland State University
% Creation Date: 2013 July 29

L = size(Sigma,1);
M = 2*L+1;

SP = (L*Sigma)^(1/2);
SP = [zeros(L,1), SP, -SP];
SP = SP + repmat(mu,1,M);

SP_x = f(SP);
P = size(SP_x,1);

mu_x = mean(SP_x,2);

Sigma_x = zeros(P);
for ix = 1:M
    t1 = SP_x(:,ix) - mu_x;
    Sigma_x = Sigma_x + t1*t1';
end
Sigma_x = Sigma_x/M;


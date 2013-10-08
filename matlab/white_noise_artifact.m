function [tau] = white_noise_artifact(u_f, u_s, a, ac, c)
% Compute the white noise artifact region on the time axis.
%
% INPUT u_f - fixed beam specified as a wave propagation vector. a B-by-3
% cartesian vector where B is the number of beams. magnitude ignored.
%
% INPUT u_s - sliding beam (see u_f). a 1-by-3 vector
%
% INPUT a - array manifold. a N-by-3 matrix of element locations
%
% INPUT ac - array manifold origin. a 1-by-3 cartesian point.
%
% INPUT c - wave speed in medium (assumed to be constant across array).
% scalar
%
% OUTPUT tau - position on time origin of delta functions comprising the
% white noise correlation artifact. a N-by-B vector.
%
% Author: John Gebbie
% Institution: Portland State University
% Creation Date: Sept 9, 2013

N = size(a,1);                                  % number of elemenbts
B = size(u_f, 1);                               % number of beams

u_f = u_f ./ repmat(sqrt(sum(u_f.^2, 2)),1,3);  % unit vectors
u_s = u_s ./ norm(u_s);
u_s = repmat(u_s, B, 1);

a = a - repmat(ac, N, 1);                       % manifold

tau = a*(u_s - u_f)' / c;

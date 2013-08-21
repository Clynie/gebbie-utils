function [R] = get_reflection_coeff(rho1, rho2, c1, c2, grz_ang)
% Compute the complex reflection coefficient for a plane boundary between
% two half space media
%
% Based on Equation (2.127) on page 95 in the following reference.
%
% F. B. Jensen, W. A. Kuperman, M. B. Porter, and H. Schmidt, Computational
% Ocean Acoustics.    New York: Springer-Verlag New York, Inc., 2000.
%
% INPUTS
%
% INPUT rho1 : density (g/cm^3) of incident medium
%
% INPUT rho2 : density (g/cm^3) of reflection medium
%
% INPUT c1 : sound speed (m/s), of incident medium
%
% INPUT c2 : sound speed (m/s), of reflection medium
%
% INPUT grz_ang : grazing angle (radians) in incident medium with respect
% to horizontal

k1 = 1./c1;
k2 = 1./c2;
kr1 = cos(grz_ang).*k1;
kz1 = sin(grz_ang).*k1;
kr2 = kr1;
kz2 = sqrt(k2.^2 - kr2.^2);

t1 = rho2.*kz1;
t2 = rho1.*kz2;
R = (t1-t2)./(t1+t2);

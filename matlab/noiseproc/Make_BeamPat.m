function [ BEAMPAT,THETA ]  = Make_BeamPat(W,FR,ZPHONES,NTHETA)

N = size(W,2);
c = 1500;
omega = 2*pi*FR; % angular FR
k = omega/c; % wavenumber
THETA = linspace(-pi/2,pi/2,NTHETA);
wpw = exp(1i*(+k*ZPHONES(:)*sin(THETA)));
BEAMPAT = nan(N,NTHETA);
for n = 1:N
    BEAMPAT(n,:) = W(:,n)' * wpw;
end

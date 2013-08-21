function V=spiral(N,R)
% V=spiral(N,R);
% N-number of points on sphere
% R- radii
% V- cartesian coordinates for
% best covering via "Generalized Spiral"

% (C) D.Tikhonov

if nargin<2,R=1;end
k=1:N;h=-1+2*(k-1)/(N-1);Tk=acos(h);
C=sqrt((8*pi/sqrt(3)));
phi=zeros(N,1);
for k=2:N-1
    x=(phi(k-1)+C/(sqrt(N)*sqrt(1-h(k)^2)));
    ix=fix(x/2/pi);
    phi(k)=x-ix*2*pi;
end
phi(N)=0;
sinTk=sin(Tk);
V=zeros(N,3);
V(:,1)=R*sinTk(:).*cos(phi(:));
V(:,2)=R*sinTk(:).*sin(phi(:));
V(:,3)=R*h(:);
end

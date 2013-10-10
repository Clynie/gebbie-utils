function [MINIMA,PATH] = gradient_descent_fcn(f, P, B, S, N)
% Simple gradient descent routine for finding local minima using function
% evals.
%
% INPUT f : function taking a M-by-1 vector and returning a scalar
%
% INPUT P : Starting point in the M-dimensional input space. a M-by-1
% vector
%
% INPUT B : bounds of P. A M-by-2 vector
%
% INPUT S : initial step size of each dimension of P. A M-by-1 vector
%
% INPUT N : number of steps
%
% OUTPUT MINIMA : position of local minima. A M-by-1 vector
%
% OUTPUT PATH : trajectory through input space
%
% Author: John Gebbie
% Creation Date: 2013 Oct 10

assert(size(P,2)==1);
M = size(P,1);
assert(size(B,1)==M);
assert(size(B,2)==2);
assert(size(S,1)==M);
assert(size(S,2)==1);

for m = 1:M
    B(m,:) = sort(B(m,:));
end

delta = S/1000;

MINIMA = P;
if nargout > 1
    PATH = nan(M,N+1);
    PATH(:,1) = P;
end

for n = 1:N
    y = nan(M,2);
    dx = delta;
    for m = 1:M
        Pd = repmat(MINIMA,1,2);
        Pd(m,:) = Pd(m,:) + [-1 1]*dx(m)/2;
        for a = 1:2
            Pd(:,a) = apply_bounds(Pd(:,a));
            y(m,a) = f(Pd(:,a));
        end
        dx(m) = diff(Pd(m,:));
    end
    dy = diff(y,[],2);
    g = dy ./ dx; % gradient
    g(dx==0) = 0; % handle dimensions that do not change
    alpha = (1-(n-1)/N).^(1/10);
    MINIMA = MINIMA - g .* S .* alpha; % minus for descent
    MINIMA = apply_bounds(MINIMA);
    if nargout > 1
        PATH(:,n+1) = MINIMA;
    end
end

    function [x] = apply_bounds(x)
        mm = x<B(:,1);
        x(mm) = B(mm,1);
        mm = x>B(:,2);
        x(mm) = B(mm,2);
    end

end


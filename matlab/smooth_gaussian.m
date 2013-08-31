function [yi] = smooth_gaussian(x,y,xi,xsigma)

%normpdf = @(x,mu,sigma) exp(-(x-mu).^2./(2.*sigma.^2));
lognormpdf = @(x,mu,sigma) -(x-mu).^2./(2.*sigma.^2);

[Xi, X] = ndgrid(xi(:), x(:));
L = lognormpdf(Xi,X,xsigma);
L = L - repmat(max(L,[],2), 1, size(L,2));
V = exp(L);
W = V ./ repmat(sum(V,2), 1, size(V,2));
yi = nan(size(xi));
yi(:) = W*y(:);

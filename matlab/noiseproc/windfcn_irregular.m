function [ Y ] = windfcn_irregular(X,WIND)

Nlowres = length(X);
Nhighres = 100*Nlowres;

Xhighres = linspace(min(X(:)),max(X(:)),Nhighres);

lowres_indices = nan(size(X));
for n = 1:Nlowres
    [~,lowres_indices(n)] = min(abs(Xhighres-X(n)));
end

Yhighres = eval(sprintf('%s(%d)',WIND,Nhighres));

Y = nan(size(X));
Y(:) = Yhighres(lowres_indices);
Y(:) = Y ./ sum(abs(Y(:))); % normalize the sum to one

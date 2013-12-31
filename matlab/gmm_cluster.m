function [Z,T,tau,mu,sigma] = gmm_cluster(x,N,iter)
% Gaussian mixture model clustering algorithm.
%
% INPUT x - input data, a K-by-D matrix consisting of K samples each of D
% dimension.
%
% INPUT N - the number of clusters
%
% INPUT iter - number of iterations
%
% OUTPUT Z - cluster assignment, a K-by-1 matrix of values in [1,N]
%
% OUTPUT T - cluster assignment strength, a K-by-N matrix of values in
% which each row sums to unity.
%
% OUTPUT tau - prior probability of each cluster, a 1-by-N matrix
%
% OUTPUT mu - means of each cluster center, a N-by-D matrix
%
% OUTPUT sigma - sigma matrices of each cluster, a D-by-D-by-N matrix
%
% Author: John Gebbie
% Institution: Metron
% Creation Date: 2013, Dec 27

warning('off','MATLAB:illConditionedMatrix')

K = size(x,1);
D = size(x,2);

sigma = repmat(diag(std(x,[],1)),[1,1,N]);
mu = repmat(std(x,[],1),N,1) .* rand(N,D) + repmat(mean(x,1),N,1);
tau = repmat(1/N,1,N);

for it = 1:iter
    f = nan(K,N);
    for n = 1:N
        f(:,n) = mvnormpdf(x,mu(n,:),sigma(:,:,n));
    end
    tauf = repmat(tau,K,1) .* f; % K-by-N
    T = tauf ./ repmat(sum(tauf,2),1,N);
    tau(:) = sum(T,1);    tau(:) = tau ./ sum(tau);
    for n = 1:N
        T2 = T(:,n) ./ sum(T(:,n));
        mu(n,:) = T2' * x;
        sigma(:,:,n) = 0;
        for k = 1:K
            x2 = x(k,:) - mu(n,:);
            sigma(:,:,n) = sigma(:,:,n) + T2(k) * (x2' * x2);
        end
    end
end

[tau,order] = sort(tau,'descend');
T = T(:,order);
mu = mu(order,:);
sigma = sigma(:,:,order);

[~,Z] = max(T,[],2);

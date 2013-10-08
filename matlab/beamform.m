function [bf] = beamform(u, dl_db, K, f, a, ac, c)
% Adaptive 3-d planewave beamforming.
%
% INPUTS
%
% INPUT u - steering vectors (beams). aka wave propagation vector. a B-by-3
% vector where B is the number of beams to form. magnitude is ignored.
%
% INPUT dl_db - diagonal loading specified in decibels relative to the mean
% power in each bin of the CSDM.
%
% INPUT K - csdm. N-by-N-by-F where N is the number of elements and F is
% the number of frequency bins
%
% INPUT f - frequencies corresponding to K. a F element vector
%
% INPUT a - array manifold. a N-by-3 matrix
%
% INPUT ac - array center. a 1-by-3 matrix
%
% INPUT c - sound speed, scalar
%
% OUTPUT bf - beamformer output. a B-by-F matrix.
%
% Author: John Gebbie
% Institution: Portland State University
% Creation Date: Sept 8, 2013

% dimensionality
B = size(u,1);                      % number of beams
N = size(a,1);                      % number of elements
F = length(f);                      % number of frequencies

% steering vectors (wave propagation directions)
u = u ./ repmat(sqrt(sum(u.^2,2)),1,3);   % B-by-3

bf = nan(B, F);                     % init output

% manifold (array element positions)
a = a - repmat(ac, N, 1);           % N-by-3

% lag distances to each element for each beam
r = a*u';                           % N-by-B

% wavenumber times -1i
minus_ik = -1i*2*pi*f/c;            % length F

parfor f1 = 1:F                        % beamform each frequency
    
    % replica vectors
    d = exp(minus_ik(f1)*r);        % N-by-B
    
    % CSDM for this frequency
    Kf = K(:, :, f1);
    
    if dl_db == Inf                 % conventional
        
        w = d;                      % N-by-B
        
    else                            % adapt
        
        % inverse CSDM
        dl = mean(diag(Kf)) * 10^(dl_db/10);    % white-noise gain
        Kwn = Kf + eye(size(Kf))*dl;            % add white noise
        
        % create adaptive steering vectors in each direction
        w = Kwn \ d;                % N-by-B
        
    end
    
    % distortionless response
    w = w ./ repmat(sum(conj(d).*w,1),N,1);     % N-by-B
    
    % beamformer kernel
    bf(:,f1) = sum((w' * Kf) .* w.', 2);
    
end

bf = real(bf);                      % toss rounding errors
bf(bf<0) = 0;

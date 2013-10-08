function [bf] = beamform_xcorr(u_f, u_s, dl_db, K, R, f, a, ac, c)
% Cross-correlate planewave beams.  Beams are distinguished as
% "fixed" and "sliding" which refer to what is conceptually
% happening to them in the time domain.
%
% INPUT u_f - wave propagation vector for "fixed" beam. a B-by-3 vector
% where B is the number of beams to form. magnitude is ignored.
%
% INPUT u_s - "sliding" beam (see u_f). a 1-by-3 vector.
%
% INPUT dl_db - diagonal loading specified in decibels relative to the mean
% power in each bin of the CSDM.
%
% INPUT K - csdm used for cross-beam correlation. N-by-N-by-F where N is
% the number of elements and F is the number of frequency bins
%
% INPUT R - csdm used for building adaptive weight vectors. N-by-N-by-F
% where N is the number of elements and F is the number of frequency bins.
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

% dimensionality
N = size(a,1);                          % number of elements
F = length(f);                          % number of frequencies
B = size(u_f,1);                        % number of beams

bf = nan(B, F);                         % init output

% steering vectors (wave propagation directions)
u_f = u_f ./ repmat(sqrt(sum(u_f.^2,2)),1,3);   % B-by-3
u_s = u_s ./ norm(u_s);                         % 1-by-3

% manifold (array element positions)
a = a - repmat(ac, N, 1);               % N-by-3

% lag distances to each element for each beam
r_f = a*u_f';                           % N-by-B
r_s = a*u_s';                           % N-by-1

% wavenumber times -1i
minus_ik = -1i*2*pi*f/c;                % length F

for f1 = 1:F                            % beamform each frequency
    
    % replica vectors
    d_f = exp(minus_ik(f1)*r_f);        % N-by-B
    d_s = exp(minus_ik(f1)*r_s);        % N-by-1
    
    % CSDM for this frequency
    Kf = K(:, :, f1);
    Rf = R(:, :, f1);
    
    if dl_db == Inf                     % conventional
        
        w_f = d_f;                      % N-by-B
        w_s = d_s;                      % N-by-1
        
    else                                % adapt
        
        % add white noise to CSDM
        dl = mean(diag(Rf)) * 10^(dl_db/10);    % white-noise gain
        Kwn = Rf + eye(size(Rf))*dl;            % add white noise
        
        % create adaptive steering vectors in each direction
        w_f = Kwn \ d_f;                % N-by-B
        w_s = Kwn \ d_s;                % N-by-1
        
    end
    
    % distortionless response
    w_f = w_f ./ repmat(sum(conj(d_f).*w_f,1),N,1);     % N-by-B
    w_s = w_s / ( d_s' * w_s );                         % N-by-1
    
    % cross-correlate beams
    bf(:,f1) = (w_f' * Kf) * w_s;
    
end


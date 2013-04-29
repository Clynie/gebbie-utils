function [PFATH,TIME,W] = Make_Pfath(CSDM,FR,FS,ZPHONES,FRMIN,FRMAX,...
    BFTYPE,ARRWIND,IFFTWIN,SHIFT,BFARGS,TRIM,TRIMFRAC,RRCDEPTHS,SOUNDSPEED)
%MAKE_PFATH         Make the passive fathometer response from noise cross
%                   spectral density.
%
%INPUT
%          CSDM ==> Cross spectral density matrix (chan, chan, freqs)
%            FR ==> Frequencies (Hz) of CSDM
%            FS ==> Sampling frequency (Hz)
%       ZPHONES ==> Vector with phones locations (starting at 0) in m.
%   FRMIN/FRMAX ==> Min/max freq (Hz)
%        BFTYPE ==> Beamforming type ('conventional', 'mvdr')
%       ARRWIND ==> i.e. 'hanning', 'taylorwin', 'boxcar' (default)
%       IFFTWIN ==> i.e. 'hanning', 'taylorwin', 'boxcar' (default)
%         SHIFT ==> Timing shift for time series misalignment
%        BFARGS ==> Arguments to beamforming algorithm
%          TRIM ==> Trim ends off of time series
%      TRIMFRAC ==> Extra trim in frac of array len (default to 0.8)
%     RRCDEPTHS ==> Depth range to normalize Rayleigh Reflection
%                   Coefficient. 2-element vector.
%    SOUNDSPEED ==> Speed of sound in water
%
%OUTPUT
%         PFATH ==> Passive fathometer derived from the noise
%          TIME ==> Time values of PFATH (sec)
%
%   Author: M. Siderius, 2/19/2009
%   Author: John Gebbie, Portland State University, 4/18/2011

if ~exist('CSDM','var'), error('must specify CSDM'); end
if ~exist('FR','var') || isempty(FR), error('must specify FR'); end
if ~exist('ZPHONES','var') || isempty(ZPHONES), error('must specify ZPHONES'); end
if ~exist('FRMIN','var') || isempty(FRMIN), FRMIN=FR(1); end
if ~exist('FRMAX','var') || isempty(FRMAX), FRMAX=FR(end); end
if ~exist('BFTYPE','var') || isempty(BFTYPE), BFTYPE='conventional'; end
if ~exist('ARRWIND','var') || isempty(ARRWIND), ARRWIND='boxcar'; end
if ~exist('IFFTWIN','var') || isempty(IFFTWIN), ARRWIND='boxcar'; end
if ~exist('SHIFT','var') || isempty(SHIFT), SHIFT=zeros(size(ZPHONES)); end
if ~exist('BFARGS','var'), BFARGS=[]; end
if ~exist('TRIM','var') || isempty(TRIM), TRIM=1; end
if ~exist('TRIMFRAC','var') || isempty(TRIMFRAC), TRIMFRAC=0.8; end
if ~exist('RRCDEPTHS','var'), RRCDEPTHS=[]; end
if ~exist('SOUNDSPEED','var') || isempty(SOUNDSPEED), SOUNDSPEED=1500; end

addpath('fftinfo');

% if number of elements specifed as time-shift of a single-element, compute
% over full array
if numel(SHIFT) == 1, SHIFT = SHIFT .* (0:length(ZPHONES)-1); end

% take subset of CSDM based on FRMIN and FRMAX
[~,idx_mnfr] = min(abs(FR - FRMIN));
[~,idx_mxfr] = min(abs(FR - FRMAX));
CSDM = CSDM(:,:,idx_mnfr:idx_mxfr);
FR = FR(idx_mnfr:idx_mxfr);

% array window function (applied spatially to array). don't assume phones
% are equally spaced. decimate a high-resolution image of the window
% function and map actual phone locations onto it.
awind = windfcn_irregular(ZPHONES(:),ARRWIND);

% pre-allocate steering vector
nchans = size(CSDM,1);
W = zeros(nchans,2,length(FR));
wc = zeros(nchans,2);

N = length(FR);
pfath_N = nan(1,N);
upbeam_N = nan(1,N);
dnbeam_N = nan(1,N);
switch (BFTYPE)
    case 'conventional'
        for n = 1:length(FR) % loop over bins
            omega = 2*pi*FR(n); % angular freq
            k = omega/SOUNDSPEED; % wavenumber
            % start with conventional steering vector
            wc(:,1) = exp(1i*(+k*ZPHONES(:)+omega*SHIFT(:)));
            wc(:,2) = exp(1i*(-k*ZPHONES(:)+omega*SHIFT(:)));
            for m = 1:2
                W(:,m,n) = awind .* wc(:,m);
                W(:,m,n) = W(:,m,n)./(wc(:,m)'*W(:,m,n)); % zero-gain
            end
            pfath_N(n) = W(:,1,n)'*CSDM(:,:,n)*W(:,2,n); % beamform and cross-correlate
            upbeam_N(n) = W(:,1,n)'*CSDM(:,:,n)*W(:,1,n);
            dnbeam_N(n) = W(:,2,n)'*CSDM(:,:,n)*W(:,2,n);
        end
    case 'mvdr'
        if isfield(BFARGS,'diag_loading')
            diag_loading = BFARGS.diag_loading;
        else
            diag_loading = 1;
        end
        csdm_real = isfield(BFARGS,'csdm_real') && BFARGS.csdm_real;
        K = zeros(nchans,nchans); % pre-allocate
        wc = zeros(nchans,2);
        for n = 1:length(FR) % loop over bins
            omega = 2*pi*FR(n); % angular freq
            k = omega/SOUNDSPEED; % wavenumber
            % invert the K matrix after adding white noise
            K(:) = CSDM(:,:,n);
            K(:) = K + eye(nchans)*mean(diag(K))*diag_loading;
            if csdm_real
                % flip array on its head to make noise field look symmetric
                % with respect to the horizontal. this forces each beam to
                % match its conjugate. this is implemented by keeping only
                % the real part of the CSDM
                K(:) = real(K);
            end
            % start with conventional steering vector
            wc(:,1) = exp(1i*(+k*ZPHONES(:)+omega*SHIFT(:)));
            wc(:,2) = exp(1i*(-k*ZPHONES(:)+omega*SHIFT(:)));
            % MVDR steering vector
            for m = 1:2
                W(:,m,n) = K\wc(:,m);
                W(:,m,n) = awind .* W(:,m,n);
                W(:,m,n) = W(:,m,n)./(wc(:,m)'*W(:,m,n)); % zero-gain
            end
            % beamform
            pfath_N(n) = W(:,1,n)'*CSDM(:,:,n)*W(:,2,n);
            upbeam_N(n) = wc(:,1)'*CSDM(:,:,n)*wc(:,1);
            dnbeam_N(n) = wc(:,2)'*CSDM(:,:,n)*wc(:,2);
        end
    otherwise
        error('unknown BFTYPE: %s', BFTYPE);
end

% apply output window on freq series prior to building full fft-series
iwind = eval(sprintf('%s(%d)',IFFTWIN,length(pfath_N))).';
pfath_N(:) = iwind .* pfath_N;

% construct a conjugate-symmetric freq-domain representation
df = abs(diff(FR(1:2)));
fftsize = round(FS/df);    % df = fs/fftsize
[fr_full,bins,~,nyqbin] = fftinfo(FS,fftsize); % bin frequencies
fr_full(nyqbin+1) = abs(fr_full(nyqbin+1)); % make nyquist positive
[~,idx_mnfr] = min(abs(fr_full - FR(1))); % min/max bin indices
[~,idx_mxfr] = min(abs(fr_full - FR(end)));
PFATH = zeros(1,fftsize);
PFATH(idx_mnfr:idx_mxfr) = pfath_N; % drop in our subset of positive freqs
pos_indices = bins>0; % make conjugate-symmetric (only positive-freq bins)
npos = sum(pos_indices);
PFATH(end:-1:end-npos+1) = conj(PFATH(pos_indices));

% bring into time-domain
PFATH(:) = real(ifft(PFATH));
TIME = 0:1/FS:(fftsize-1)/FS;

if TRIM
    % mask out garbage spanning 2X length of array at start (and leaking
    % into back of time series
    len = 2*abs(max(ZPHONES(:))-min(ZPHONES(:)));
    len_extra_buffer = TRIMFRAC * len;
    mask_a = min(fftsize,ceil((len+len_extra_buffer)/SOUNDSPEED*FS));
    mask_b = fftsize-min(fftsize,ceil(len_extra_buffer/SOUNDSPEED*FS));
    PFATH = PFATH(mask_a:mask_b);
    TIME = TIME(mask_a:mask_b);
end

% compute the relative power of up and down beams
if isempty(RRCDEPTHS) || length(RRCDEPTHS) ~= 2
    rrca = 1;
    rrcb = length(PFATH);
else
    depth = TIME*c/2;
    [~,rrca] = min(abs(depth-RRCDEPTHS(1)));
    [~,rrcb] = min(abs(depth-RRCDEPTHS(2)));
end
reflection_energy_ratio = sum(abs(upbeam_N).^2,2) / sum(abs(dnbeam_N).^2,2);
PFATH(:) = sqrt(reflection_energy_ratio / sum(abs(PFATH(rrca:rrcb)).^2)) * PFATH;

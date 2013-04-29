function [RL,RL_THETA,BF,BF_THETA] = Make_RL(CSDM,FR,ZPHONES,NANG,...
    BFTYPE,ARRWIND,SHIFT,BFARGS)
%MAKE_RL            Make the reflection loss from noise cross spectral
%                   density.
%
%INPUT
%          CSDM ==> Cross spectral density matrix (chan, chan, freqs)
%            FR ==> Frequencies (Hz)
%       ZPHONES ==> Vector with phones locations (starting at 0) in m.
%          NANG ==> Number of RL angles (defaults to 100)
%        BFTYPE ==> Beamforming type ('conventional', 'mvdr')
%       ARRWIND ==> i.e. 'hanning', 'taylorwin', 'boxcar' (default)
%         SHIFT ==> Timing shift for time series misalignment
%        BFARGS ==> Arguments to beamforming algorithm
%
%OUTPUT
%            RL ==> Reflection loss derived from the noise (angs, freqs)
%      RL_THETA ==> Angles for reflection loss (degrees)
%            BF ==> Beamformer output (angs, freqs)
%      BF_THETA ==> Angles for beamformer output (degrees)
%
%NOTES
%                   This function preserves the signal energy of the CSDM.
%                   The BF output signal at *each* angle should be on the
%                   order of a *single* time-series input used to create
%                   the CSDM in the first place.
%
%   Author: M. Siderius, HLS Research Inc., 9/6/2007
%   Author: John Gebbie, Portland State University, 4/7/2011

if ~exist('CSDM','var'), error('must specify CSDM'); end
if ~exist('FR','var') || isempty(FR), error('must specify FR'); end
if ~exist('ZPHONES','var') || isempty(ZPHONES), error('must specify ZPHONES'); end
if ~exist('NANG','var') || isempty(NANG), NANG=100; end
if ~exist('BFTYPE','var') || isempty(BFTYPE), BFTYPE='conventional'; end
if ~exist('ARRWIND','var') || isempty(ARRWIND), ARRWIND='boxcar'; end
if ~exist('SHIFT','var') || isempty(SHIFT), SHIFT=zeros(size(ZPHONES)); end
if ~exist('BFARGS','var'), BFARGS=[]; end

% if number of elements specifed as time-shift of a single-element, compute
% over full array
if numel(SHIFT) == 1, SHIFT = SHIFT .* (0:length(ZPHONES)-1); end

% allocate return matrices and angle vectors
BF_THETA = linspace(90,-90,NANG*2+1);
RL_THETA = linspace(0,90,NANG);
BF = zeros(length(BF_THETA),length(FR));
RL = zeros(length(FR),length(RL_THETA));

nbfang = length(BF_THETA); % number of beamforming angles
nchans = length(ZPHONES); % number of channels

% array window function (applied spatially to array). don't assume phones
% are equally spaced. decimate a high-resolution image of the window
% function and map actual phone locations onto it.
wind = windfcn_irregular(ZPHONES,ARRWIND).';

% speed of sound in water
c = 1500;

% pre-compute some values used inside beamforming loops
dsintheta = ZPHONES(:)*sin(BF_THETA*pi/180); % (phones, theta)

% pre-allocate steering vector
w = zeros(nchans,1);
wConv = zeros(nchans,1);

switch (BFTYPE)
    case 'conventional'
        for n = 1:length(FR) % loop over bins
            omega = 2*pi*FR(n); % angular freq
            k = omega/c; % wavenumber
            for m = 1:nbfang % loop over beams
                % start with conventional steering vector
                wConv(:) = exp(1i*(k*dsintheta(:,m)+omega*SHIFT(:)));
                w(:) = wind(:) .* wConv;
                w(:) = w./(wConv'*w); % normalize the sum to one
                BF(m,n) = w'*CSDM(:,:,n)*w; % beamform
            end
        end
    case { 'mvdr' }
        Kinv = zeros(nchans,nchans); % pre-allocate
        if isfield(BFARGS,'diag_loading')
            diag_loading = BFARGS.diag_loading;
        else
            diag_loading = 1;
        end
        if isfield(BFARGS,'csdm_symmetric')
            csdm_symmetric = BFARGS.csdm_symmetric;
        else
            csdm_symmetric = 0;
        end
        if isfield(BFARGS,'csdm_real')
            csdm_real = BFARGS.csdm_real;
        else
            csdm_real = 0;
        end
        for n = 1:length(FR) % loop over bins
            omega = 2*pi*FR(n); % angular freq
            k = omega/c; % wavenumber
            % invert the K matrix after adding white noise
            Kinv(:) = CSDM(:,:,n);
            if csdm_symmetric
                % flip array on its head to make noise field look symmetric
                % with respect to the horizontal. this forces each beam to
                % match its conjugate
                Kinv(:) = Kinv + rot90(flipud(Kinv));
            end
            if csdm_real
                % keep only the real part of the CSDM
                Kinv(:) = real(Kinv);
            end
            Kinv(:) = Kinv + eye(nchans)*mean(diag(Kinv))*diag_loading;
            Kinv(:) = inv(Kinv);
            for m = 1:nbfang % loop over beams
                % start with conventional steering vector
                wConv(:) = exp(1i*(k*dsintheta(:,m)+omega*SHIFT(:)));
                w(:) = wind(:) .* wConv;
                w(:) = Kinv*w/(w'*Kinv*w); % MVDR steering vector
                w(:) = w./(wConv'*w); % normalize the sum to one
                BF(m,n) = w'*CSDM(:,:,n)*w; % beamform
            end
        end
    otherwise
        error('unknown BFTYPE: %s', BFTYPE);
end

BF(:) = abs(BF); % toss out imaginary parts due to rounding
% compute RL by dividing the negative-angle beams (sound coming from
% seabed) by the positive-angle beams (sound coming from surface). this
% *should* yield a value that is less than one since seabed absorbs energy.
RL(:) = rot90( ( BF(NANG+2:end,:) + eps ) ./ ( BF(NANG:-1:1,:) + eps ) );

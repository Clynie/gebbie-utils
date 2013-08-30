function [CSDM,FR,NSAMP,FR_FULL] = Make_CSDMatrix(TSERIES,FS,FFTSIZE,...
    FRMIN,FRMAX,WIND,OVERLAP,ZEROPAD,PREWHITEN)
%MAKE_CSDMATRIX         Cross Spectral Density Matrix
%
%INPUT
%           TSERIES ==> time series data, each channel in a separate row
%                FS ==> sampling frequency in Hz
%           FFTSIZE ==> number of samples in each time interval
%       FRMIN/FRMAX ==> minimum and maximum frequencies for which to
%                       compute CSDM. If not specified, or specified as [],
%                       defaults to 0 and FS/2, respectively.
%              WIND ==> window function to apply to each time interval. If
%                       not specified, or specified as [], defaults to
%                       'hanning'. To turn off windowing, use 'boxcar'.
%           OVERLAP ==> fractional overlap of each time interval. If not
%                       specified, or specified as [], defaults to .5.
%           ZEROPAD ==> zero-pad the end of each buffer input to the FFT by
%                       padding the input buffer. If is a boolean set to
%                       true, pads by twice the size. If a positive integer
%                       greater than one, pads to that exact number of
%                       samples.
%         PREWHITEN ==> pre-whiten the FFT output by setting the magnitudes
%                       of each frequency bin to one while retaining the
%                       complex phase.  Note that this will destroy the
%                       relative and absolute power levels of the input
%                       signals.
%
%OUTPUT
%              CSDM ==> N-by-N-by-M CSDM where N is the number of channels
%                       and M is the number of bins between FRMIN and
%                       FRMAX.
%                FR ==> the frequencies associated with each bin
%             NSAMP ==> the number of samples averaged into the CSDM
%           FR_FULL ==> the full set of frequencies for the expanded window
%
%   Author: M. Siderius, HLS Research Inc., 9/6/2007
%   Author: John Gebbie, Portland State University, 4/7/2011

if ~exist('TSERIES','var'), error('must specify TSERIES'); end
if ~exist('FS','var') || isempty(FS), error('must specify FS'); end
if ~exist('FFTSIZE','var') || isempty(FFTSIZE), error('must specify FFTSIZE'); end
if ~exist('FRMIN','var') || isempty(FRMIN), FRMIN=0; end
if ~exist('FRMAX','var') || isempty(FRMAX), FRMAX=FS/2; end
if ~exist('WIND','var'), WIND='hanning'; end
if ~exist('OVERLAP','var'), OVERLAP=.5; end
if ~exist('ZEROPAD','var') || isempty(ZEROPAD), ZEROPAD=false; end
if ~exist('PREWHITEN','var') || isempty(PREWHITEN), PREWHITEN=false; end

if isnumeric(ZEROPAD) && ZEROPAD > 1
    FFTSIZE_IN = ZEROPAD;
elseif ZEROPAD
    FFTSIZE_IN = 2*FFTSIZE;
else
    FFTSIZE_IN = FFTSIZE;
end

% the frequencies of each bin of the DFT output
FR_FULL = fftinfo(FS,FFTSIZE_IN,true);

% bin indices of min and max frequences to compute CSDM for
[~,frmin_ix] = min(abs(FR_FULL-FRMIN));
[~,frmax_ix] = min(abs(FR_FULL-FRMAX));
FR = FR_FULL(frmin_ix:frmax_ix);

nchans = size(TSERIES,1); % number of channels
CSDM = zeros(nchans,nchans,length(FR)); % CSDM starts at zero
fftin = zeros(nchans,FFTSIZE_IN); % tmp buffer to hold fft input
fftout = nan(nchans,FFTSIZE_IN); % tmp buffer to hold fft output
p = zeros(nchans,1); % tmp buffer to hold vector of pressure values

% window applied to time series interval prior to FFT (prevents ringing
% caused by abrupt edges at edge of intervals)
timewindfcn = eval(sprintf('%s(%d)',WIND,FFTSIZE)).';
% normalize power of window function to unity. this preserves the total
% energy of the time series to which it is multiplied.
timewindfcn = timewindfcn.*sqrt(FFTSIZE/(sum(abs(timewindfcn).^2)));
timewindfcn = repmat(timewindfcn,nchans,1);

% overlap time windows by 50%
intvl_stride = max(1,FFTSIZE*(1-OVERLAP));
intvl_start_indices = round(1:intvl_stride:size(TSERIES,2)-FFTSIZE+1);
for m = intvl_start_indices
    fftin(:,1:FFTSIZE) = TSERIES(:,m:m+FFTSIZE-1);
    fftin(:,1:FFTSIZE) = fftin(:,1:FFTSIZE)-repmat(mean(fftin(:,1:FFTSIZE),2),1,FFTSIZE); % zero-mean
    fftin(:,1:FFTSIZE) = fftin(:,1:FFTSIZE).*timewindfcn; % apply time window
    % apply time window, and FFT each channel separately
    fftout(:) = fft(fftin,[],2);
    if PREWHITEN
        fftout = fftout ./ abs(fftout);
    end
    % calcualte CSDM for each bin separately
    for n = 1:length(FR)
        p(:) = fftout(:,n+frmin_ix-1); % data values for each array element
        CSDM(:,:,n) = CSDM(:,:,n) + p*p'; % hermitian outter prod, accum.
    end
end

% the number of samples that were accumulated into the CSDM
nsamp_avg = length(intvl_start_indices)*FFTSIZE;
% the number of samples of wall-clock time (non-repeated samples) that were
% accumulated into the CSDM
NSAMP = intvl_start_indices(end)+FFTSIZE-1;

% scale the CSDM so that the energy is proportional to the amount of
% wall-clock time that was actually averaged
scale = NSAMP/nsamp_avg;
% scale the CSDM by the number of individual CSDMs added together
scale = scale/length(intvl_start_indices);
% scale the CSDM so that each bin is in units of power/hz
scale = scale/(FS*FFTSIZE);
% apply the scaling
CSDM = CSDM .* scale;

classdef Freq < handle
    % Utility class for dealing with single dimension contiguous band
    % frequency-domain data.
    %
    % Author: John Gebbie
    % Creation Date: 2013 Apr 08
    
    properties (SetAccess = private)
        fs      % sample rate (Hz)
        N       % fftsize (total number of bins)
        M       % number of bins in partial spectrum (fr)
        band    % a 2-element vector of min/max frequencies (inclusive)
        
        full    % a N-by-1 vector of all frequencies
        mask    % a N-by-1 boolean vector mask of 'full'
        fr      % a M-by-1 vector of partial frequencies, in which M < N
    end
    
    
    methods % public
        
        function [f] = synthFull(o, part)
            % Synthesize full spectrum
            %
            % INPUT part : input partial spectrum (in columns)
            %
            % OUTPUT f : full spectrum
            
            % ensure the partial spectrum has the proper dimensions by
            % compressing to a M-by-Nch matrix
            S = size(part);
            assert(o.M == length(o.fr));
            assert(o.M == S(1));
            part = reshape(part, o.M, []);
            
            % synthesize full spectrum
            f = fftsynth(part, o.fr(1), o.full);
            
            % restore higher dimensionality
            S(1) = o.N;
            f = reshape(f, S);
        end
        
        function [t] = timeAxis(o)
            % Time axis for this snapshot size
            %
            % OUTPUT t : N-by-1 vector of time offsets, starting at t=0
            
            t = linspace(0, (o.N-1)/o.fs, o.N).';
        end
        
        function [y, t] = synthTime(o, part, wrap, env, wfn)
            % Synthesize partial spectrum into time-series
            %
            % INPUT part : partial spectrum (in columns)
            %
            % INPUT wrap : wrap the time-series so t=zero is in the center
            % and the last half is wrapped to negative time.
            %
            % INPUT env : compute envelope
            %
            % INPUT wfn : apply hanning window to frequency domain
            % (boolean)
            %
            % OUTPUT y : the output time series (in columns), having
            % dimension N-by-C in which C is the number of channels.
            %
            % OUTPUT t : time time axis a N-by-1 column vector
            
            % process input args
            if nargin < 3 || isempty(wrap)
                wrap = false;
            end
            if nargin < 4 || isempty(env)
                env = false;
            end
            if nargin < 5 || isempty(wfn)
                wfn = false;
            end
            
            % ensure the partial spectrum has the proper dimensions by
            % compressing to a M-by-Nch matrix
            S = size(part);
            assert(o.M == length(o.fr));
            assert(o.M == S(1));
            part = reshape(part, o.M, []);
            Nch = S(2);
            
            % apply window function
            if wfn
                part = part .* repmat(hanning(o.M), 1, Nch);
            end
            
            % synthesize full spectrum
            f = o.synth(part);
            
            % IFFT
            y = real(ifft(f));
            t = o.timeAxis();
            
            % wrap the time axis
            if wrap
                [t, y] = wrap_time(t, y);
            end
            
            % envelope the signal
            if env
                y = abs(hilbert(y));
            end
            
            % restore higher dimensionality
            S(1) = o.N;
            y = reshape(y, S);
        end
        
        function [bin_fr, full_idx, fr_idx] = get_nearest_bin(o, in_fr)
            % Nearest bin to the specified frequency
            %
            % INPUT in_fr - input frequency (Hz)
            %
            % OUTPUT bin_fr : the frequency of the bin
            %
            % OUTPUT full_idx : bin index.  this is the index into the
            % 'full' spectrum corresponding that is closest to the input
            % frequency
            %
            % OUTPUT fr_idx : same as full_idx, but index in the 'fr'
            % spectrum array.
            
            [~, fr_idx] = min(abs(in_fr - o.fr));
            [~, full_idx] = min(abs(in_fr - o.full));
            assert(o.fr(fr_idx) == o.full(full_idx)); % sanity check
            bin_fr = o.full(full_idx);
        end
        
        function [freq] = newScaledFreq(o, factor)
            % Create a new Freq object that scales the sample rate an
            % integer amount.  This facilitates fast upsampling since the
            % same partial spectrum can be used to generate a time series
            % with the returned Freq object.
            %
            % INPUT factor : positive integer scale factor
            %
            % OUTPUT freq : a new Freq object.
            
            assert(factor == floor(factor) && factor >= 1, ...
                'factor must be a positive integer');
            freq = Freq.newBasic(o.fs*factor, o.N*factor, o.fr([1 end]));
        end
        
        function [s, t] = spectrogram(o, y, wind, overlap)
            % Create a spectrogram from a long time series using. Segments
            % the time series into length N snapshots and computes FFT of
            % each one. Returns result in a matrix easily displayed as a
            % spectrogram. Optionally plots a spectrogram if no outputs are
            % requested.  To plot a spectrogram manually, one would write:
            %
            %       imagesc(o.fr, t, 20*log10(abs(s)));
            %
            % INPUT y : the time series. Should be single dimensional
            %
            % INPUT wind : window function applied to each snapshot. Should
            % be 1-by-N in which N is defined in this class. Will be
            % rescaled to prevent adding power to the input time series.
            %
            % INPUT overlap : overlap of snapshots. Must be between 0
            % (inclusive) and 1 (exclusive). A value closer to 0 means
            % greater overlap, and thus more snapshots. Recommend using a
            % value greater than 0.1. Defaults to 0.5.
            %
            % OUTPUT s : P-by-M in which P is the number of snapshots and M
            % is defined in this class as the number of bins. Returns
            % complex outputs of FFT, with the squred-magnitude in units of
            % power per Hz.
            %
            % OUTPUT t : P-by-1 of time offsets of each snapshot.
            
            narginchk(2,4);
            if nargin < 3 || isempty(wind)
                wind = ones(o.N,1);
            end
            if nargin < 4
                overlap = 0.5;
            end
            
            assert(overlap >= 0 && overlap < 1, ...
                'overlap must be between 0 (inclusive) and 1 (exclusive)');
            assert(numel(wind) == o.N, ...
                'wind must have exactly N (%d) elements', o.N);
            
            snap_idx = floor(1 : o.N*(1-overlap) : numel(y)-o.N+1).';
            t = (snap_idx-1)/o.fs;  % time axis

            P = length(snap_idx);
            idx = repmat(snap_idx, 1, o.N) + repmat((0:o.N-1), P, 1);
            s = nan(size(idx));
            s(:) = y(idx(:));
            wind = sqrt(o.N/sum(abs(wind).^2)).*wind;
            s = s .* repmat(wind(:)', [P 1]);
            s = fft(s, [], 2);
            s = s(:, o.mask) ./ sqrt(o.fs*o.N) .* 2;
            
            if nargout == 0
                imagesc(o.fr, t, 20*log10(abs(s)));
            end
        end
        
    end
    
    
    methods (Access = private)
        
        function [o] = Freq(fs, N, band)
            % Constructor
            %
            % INPUT fs : sample rate (Hz)
            %
            % INPUT N : number of FFT points
            %
            % INPUT band : a 2-element vector of positive frequencies
            % (optional -- defaults to all positive frequencies)
            
            % process input args
            if nargin < 3 || isempty(band)
                band = [eps Inf];
            end
            assert(all(band >= 0));
            
            % save basic parameters
            o.fs = fs;
            o.N = N;
            o.band = band;
            
            % compute the full spectrum
            o.full = Freq.fftinfo(o.fs, o.N).';
            % conpute mask
            o.mask = o.band(1) <= o.full & o.full <= o.band(2);
            
            % handle case when no frequency bins fall within the band
            % limits since this is probably just a single frequency
            if ~any(o.mask)
                [~, idx] = min(abs(o.full - mean(band)));
                o.mask(idx) = true;
            end
            
            % compute partial spectrum
            o.fr = o.full(o.mask);
            % count frequency bins
            o.M = length(o.fr);
        end
        
    end
    
    
    methods (Static)
        
        function [freq] = newBasic(fs, N, band)
            % Create a new Freq object using a fixed frequency resolution
            % (i.e. width of each frequency bin).
            %
            % INPUT fs : sample rate (Hz)
            %
            % INPUT N : number of FFT points
            %
            % INPUT band : a 2-element vector of positive frequencies
            % (optional -- defaults to all positive frequencies)
            
            % process input args
            if nargin < 3
                band = [];
            end
            
            freq = Freq(fs, N, band);
        end
        
        function [freq] = newBinRes(fs, res, band)
            % Create a new Freq object using a fixed frequency resolution
            % (i.e. width of each frequency bin).
            %
            % INPUT fs : sample rate (Hz)
            %
            % INPUT res : width of each frequency bin (Hz)
            %
            % INPUT band : a 2-element vector of positive frequencies
            % (optional -- defaults to all positive frequencies)
            
            % process input args
            if nargin < 3
                band = [];
            end
            
            N1 = max(round(fs / res), 1);
            freq = Freq(fs, N1, band);
        end
        
        function [freq] = newByTime(fs, T, band)
            % Create a new Freq object using a fixed frequency resolution
            % (i.e. width of each frequency bin).
            %
            % INPUT fs : sample rate (Hz)
            %
            % INPUT T : duration of time window (s)
            %
            % INPUT band : a 2-element vector of positive frequencies
            % (optional -- defaults to all positive frequencies)
            
            % process input args
            if nargin < 3
                band = [];
            end
            
            N1 = max(round(T * fs), 1);
            freq = Freq(fs, N1, band);
        end
        
        function [freqs] = fftinfo(fs, N, nyqpos)
            % Bin frequencies output by fft
            %
            % INPUT fs : Sampling frequency
            %
            % INPUT N : Number of samples output by fft. This is specified
            % by the 2nd paramter to 'fft'.
            %
            % INPUT nyqpos : <optional> Boolean. Make the nyquist frequency
            % positive for even numbered N's. This frequency is negative by
            % default.
            %
            % OUTPUT freqs : Frequency of each bin (1-by-N)
            
            if nargin < 3
                nyqpos = false;
            end
            
            binwidth = fs/N;
            bins = (mod(((0:N-1)+floor(N/2)), N)-floor(N/2));
            nyq = abs(bins(floor(N/2)+1));
            if nyqpos && mod(N,2) == 0
                bins(nyq+1) = abs(bins(nyq+1));
            end
            freqs = binwidth.*bins;
        end
        
        
    end
end

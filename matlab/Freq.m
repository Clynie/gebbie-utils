classdef Freq < handle
    % Utility class for handling single dimension frequency-domain data.
    %
    % Author: John Gebbie
    % Creation Date: 2013-04-08
    
    properties (SetAccess = private)
        fs      % sample rate (Hz)
        N       % fftsize
        band    % a 2-element vector of min/max frequencies (inclusive)
        
        full    % a N-by-1 vector of all frequencies
        mask    % a N-by-1 boolean vector mask of 'full'
        fr      % a M-by-1 vector of partial frequencies, in which M < N
    end
    
    
    methods % public
        
        function [f] = synth(o, part)
            % synthesize full spectrum
            %
            % INPUT part        : input partial spectrum (in columns)
            %
            % OUTPUT f          : full spectrum
            
            % ensure the partial spectrum has the proper dimensions by
            % compressing to a M-by-Nch matrix
            S = size(part);
            M = size(part, 1);
            assert(M == length(o.fr));
            part = reshape(part, M, []);
            
            % synthesize full spectrum
            f = fftsynth(part, o.fr(1), o.full);
            
            % restore higher dimensionality
            S(1) = o.N;
            f = reshape(f, S);
        end
        
        function [t] = timeAxis(o)
            % get time axis for this snapshot size
            %
            % OUTPUT t          : N-by-1 vector of time offsets, starting
            % at t=0
            
            t = linspace(0, (o.N-1)/o.fs, o.N).';
        end
        
        function [y, t] = synthTime(o, part, wrap, env, wfn)
            % synthesize partial spectrum into time-series
            %
            % INPUT part        : partial spectrum (in columns)
            %
            % INPUT wrap        : wrap the time-series so t=zero is in the
            % center and the last half is wrapped to negative time.
            %
            % INPUT env         : compute envelope
            %
            % INPUT wfn         : apply hanning window to frequency domain
            % (boolean)
            %
            % OUTPUT y          : the output time series (in columns),
            % having dimension N-by-C in which C is the number of channels.
            %
            % OUTPUT t          : time time axis a N-by-1 column vector
            
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
            M = size(part, 1);
            assert(M == length(o.fr));
            part = reshape(part, M, []);
            Nch = size(part, 2);
            
            % apply window function
            if wfn
                part = part .* repmat(hanning(M), 1, Nch);
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
            % Returns the nearest bin to the specified frequency.
            %
            % INPUTS
            %
            %   in_fr           : input frequency (Hz)
            %
            % OUTPUTS
            %
            %   bin_fr          : the frequency of the bin
            %
            %   full_idx        : bin index.  this is the index into the
            %   'full' spectrum corresponding that is closest to the input
            %   frequency
            %
            %   fr_idx          : same as full_idx, but index in the 'fr'
            %   spectrum array.
            
            [~, fr_idx] = min(abs(in_fr - o.fr));
            [~, full_idx] = min(abs(in_fr - o.full));
            assert(o.fr(fr_idx) == o.full(full_idx)); % sanity check
            bin_fr = o.full(full_idx);
        end
        
    end
    
    
    methods (Access = private)
        
        function [o] = Freq(fs, N, band)
            % Constructor
            %
            % INPUT fs      : sample rate (Hz)
            %
            % INPUT N       : number of FFT points
            %
            % INPUT band    : a 2-element vector of positive frequencies
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
            o.full = fftinfo(o.fs, o.N).';
            % conpute mask
            o.mask = o.band(1) <= o.full & o.full <= o.band(2);
            % compute partial spectrum
            o.fr = o.full(o.mask);
        end
        
    end
    
    
    methods (Static)
        
        function [freq] = newBasic(fs, N, band)
            % Create a new Freq object using a fixed frequency resolution
            % (i.e. width of each frequency bin).
            %
            % INPUT fs      : sample rate (Hz)
            %
            % INPUT N       : number of FFT points
            %
            % INPUT band    : a 2-element vector of positive frequencies
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
            % INPUT fs      : sample rate (Hz)
            %
            % INPUT res     : width of each frequency bin (Hz)
            %
            % INPUT band    : a 2-element vector of positive frequencies
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
            % INPUT fs      : sample rate (Hz)
            %
            % INPUT T       : duration of time window (s)
            %
            % INPUT band    : a 2-element vector of positive frequencies
            % (optional -- defaults to all positive frequencies)
            
            % process input args
            if nargin < 3
                band = [];
            end
            
            N1 = max(round(T * fs), 1);
            freq = Freq(fs, N1, band);
        end
        
    end
end

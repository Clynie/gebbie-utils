function [ freqs,bins,binwidth,nyq ] = fftinfo(sampfreq,fftsize,nyqpos)
% FFTINFO Produce useful information about a fft window.
%   INPUT
%       sampfreq        => sampling frequency
%       fftsize         => size of vector *output* by fft. this is
%                               specified by the 2nd paramter to 'fft'.
%       nyqpos          => make the nyquist frequency positive for even
%                               numbered fftsizes. by default it is
%                               negative.
%   OUTPUT
%       freqs           => the frequency of each bin (1-by-fftsize)
%       bins            => the bin numbers (1-by-fftsize)
%       binwidth        => the width of each bin (scalar)
%       nyq             => the nyquest frequency (scalar)
%
% Author: John Gebbie

binwidth = sampfreq/fftsize;
bins=(mod(((0:fftsize-1)+floor(fftsize/2)), fftsize)-floor(fftsize/2));
nyq=abs(bins(floor(fftsize/2)+1));
if exist('nyqpos','var') && nyqpos && mod(fftsize,2)==0
    bins(nyq+1) = abs(bins(nyq+1));
end
freqs=binwidth.*bins;

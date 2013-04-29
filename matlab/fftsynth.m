function [ spect_full ] = fftsynth(spect_partial,freq_low,freqs_full,dim)
% FFTSYNTH Expand signal specified with only positive frequencies to
% include the full set of positive/negative bins. The output can be then
% passed to ifft to produce a real signal.
%   INPUT
%       spect_partial  => input signals (selected contiguous bins)
%       freq_low        => lowest frequency in spect_partial. if specified
%                           as a vector, only the first element is used.
%       freqs_full      => full set of frequency bins to output (see
%                           fftinfo)
%       dim             => dimension along which signals are laid out. if
%                           signals are in columns, dim is 1. default is 1.
%   OUTPUT
%       spect_full     => signal including full set of frequency bins
%
% Author: John Gebbie

if ~exist('dim','var'), dim=1; end

if dim==2, spect_partial = spect_partial.'; end

num_chans = size(spect_partial,2);
% find the first bin spect_partial will fit into
frmin_idx = find(freq_low(1)==freqs_full,1,'first');
assert(~isempty(frmin_idx),'freq_low does not match any frequencies in freqs_full');
% zero-pad frequencies leading up to first bin
spect_full = [ zeros(frmin_idx-1,num_chans) ; spect_partial ];
% find the nyquist bin. if freqs_full is even, nyquist bin is negative, so
% use abs
[~,frnyq_idx] = max(abs(freqs_full));
assert(size(spect_full,1) <= frnyq_idx,'supplied spect_partial has more bins than the positive frequencies in freqs_full');
% zero-pad up to the nyquist bin (the first nyquist bin in the case
% freqs_full is of odd length
spect_full = [ spect_full ; zeros(frnyq_idx-size(spect_full,1),num_chans) ];
% all that should remain are the negative frequency bins
num_negative_bins = length(freqs_full) - size(spect_full,1);
% bin 0 is DC, so start at idx 2. for real signals, the negative bins are
% mirrored over the middle (nyquist) bin and conjugated.
spect_full = [ spect_full ; conj(spect_full(num_negative_bins+1:-1:2,:)) ];

if dim==2, spect_full = spect_full.'; end

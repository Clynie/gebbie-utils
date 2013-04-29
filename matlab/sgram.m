function [spect, fr_full, fr_mask, time_axis, freq_axis, slot_indices] = ...
    sgram(tseries, fs, fr_res, fr_band, wind, overlap_frac)

if ~exist('fr_band', 'var'), fr_band = []; end
if ~exist('overlap_frac', 'var') || isempty(overlap_frac), ...
        overlap_frac = []; end

fftsize = max(1, round(fs/fr_res));
fr_full = fftinfo(fs, fftsize, true);
if ~isempty(fr_band)
    fr_mask = fr_band(1) <= fr_full & fr_full <= fr_band(2);
else
    fr_mask = 0 < fr_full;
end

slot_indices_start = (1:ceil(fftsize*(1-overlap_frac)):size(tseries, 1)-fftsize+1).';
nslots = length(slot_indices_start);
slot_indices = repmat(slot_indices_start, 1, fftsize) + ...
    repmat((0:fftsize-1), nslots, 1);
spect = nan(size(slot_indices));
spect(:) = tseries(slot_indices(:));
wfn = eval(sprintf('%s(%d)',wind,fftsize)).';
wfn = sqrt(fftsize/sum(abs(wfn).^2)).*wfn;
spect = spect .* repmat(wfn, [nslots 1]);
spect = fft(spect, [], 2);
spect = spect(:, fr_mask) ./ sqrt(fs*fftsize) .* 2;

freq_axis = fr_full(fr_mask);
time_axis = (slot_indices_start-1)/fs;

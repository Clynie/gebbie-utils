function [y] = transfer_fcn_to_sim(tf, freq, out_smp, ...
    source_db, noise_db)
% Synthesize acoustic noise data using transfer function.
%
% INPUT tf :: transfer function Nb-by-Ne in which Nb is the number of
% bins and Ne is the number of elements
%
% INPUT freq :: Freq object corresponding to the bins in tf. freq.M == Nb
%
% INPUT nsmp :: number of output samples
%
% INPUT source_db, noise_db :: source level and noise level (this is not
% SNR since it does not incorporpote transmission loss)
%
% OUTPUT y :: output time series.  Ns-by-Ne matrix in which Ns == out_smp
%
% Author: John Gebbie
% Institution: Portland State University
% Creation Date: 2013-09-06

assert(size(tf, 1) == freq.M);
Ne = size(tf, 2);

S = repmat(exp(1i*2*pi*rand(freq.M, 1)), 1, Ne);
N = exp(1i*2*pi*rand(freq.M, Ne));

f = S.*tf.*10^(source_db/20) + N.*10^(noise_db/20);

y = freq.synthTime(f, false, false, false);
y = y(1:out_smp, :);

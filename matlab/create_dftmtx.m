function [dftm] = create_dftmtx(fftsize,fr,fs)
% Create a matrix to compute a discrete Fourier transform (DFT)
%
% INPUT fftsize : FFT size (number of samples in time series)
%
% INPUT fr : column vector of N frequencies to compute DFT at
%
% INPUT fs : sampling frequency of timeseries
%
% OUTPUT dftm : a N-by-fftsize matrix. To compute DFT, left multiply a
% timeseries matrix by this matrix. The timeseries matix has its series
% values in columns with each column a separate series.

dftm=exp(-2*pi*1i*fr./fs*(0:fftsize-1));

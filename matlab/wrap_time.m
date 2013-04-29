function [time,signal] = wrap_time(time,signal)
N = length(time);
S = ceil(N/2);
T = ( time(end)-time(1) ) / (N-1) * N;
time(S:N) = time(S:N)-T;
timemask = [S:N 1:S-1];
time = time(timemask);
signal = signal(timemask,:);

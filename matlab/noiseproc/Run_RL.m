%% init
clear;
clc;


%% read data and build cross spectral density matrices

dsize = 2^7;
dirname='J:\disk1_backup\MFA2003\boundary03-mfa-float-an10';
[name,nfiles] = get_fnames(dirname);        %get names of files (specific for NURC files) needs to be changed depending on data
n=220;
mxsamp=90*12000;
nsamp=0;
zphones = (0:.18:5.58)'; %array geometry
chans = length(zphones);

%loop over some number of files (say 1 minute of data) to form a averaged
%cross spectral density matrix
while nsamp<mxsamp
    fprintf('processing %s (%3d%% done)\n',name(n,:),floor(nsamp/mxsamp*100));
    %read in data
    out = raw([dirname filesep() name(n,:)]);     %read in the binary data (NURC format)
    fs = out.samp_f;
    %make cross spectral density matrix
    [covtmp,fr,nsamp1] = Make_CSDMatrix_RL(out.t_s,fs,dsize);
    nsamp=nsamp+nsamp1;
    %average csdm's over files
    if(~exist('cov','var'))
        cov =  covtmp;
    else
        cov = cov + covtmp;
    end
    n=n+1;
end


%% beam form

%beamform and make reflection loss
shift = 2.6e-6;
[Noise_RL,theta,bf_out,phi] =  Make_RL(cov,fr,zphones,1,'H',shift);


%% make plot

close all;

%image beamformer output
figure;
pcolor(fr,phi,10*log10(bf_out));shading('flat');axis([100,6000,-90,90,]);
title('Beam-former Output (units unknown)');
xlabel('Frequency (Hz)');
ylabel('Beam angle (deg)');
colorbar;

%image RL
figure;
pcolor(theta,fr,Noise_RL');shading('flat');axis([0,90,200,4000]);caxis([0,15])
title('Reflection Loss (dB)');
xlabel('Grazing angle (deg)');
ylabel('Frequency (Hz)');
colorbar;

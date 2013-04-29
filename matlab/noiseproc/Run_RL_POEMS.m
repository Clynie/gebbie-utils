%% init
clear;
clc;
addpath('..\poems\poems_a_parser');
close all;
figure; set(gcf,'WindowStyle','docked');

idirpath='J:\vla\2010_10_29_13_00_00_deployment2';
ifilename='PF05_2010_10_30_07_*_01.hti';
ifilestructs=dir(fullfile(idirpath,ifilename));

for ixfile=1:length(ifilestructs)
    
    %% read data
    
    fprintf('reading %s...\n',ifilestructs(ixfile).name);
    ifilepath=fullfile(idirpath,ifilestructs(ixfile).name);
    tseries = poems_a_hti(ifilepath,1,[],[],'double',0,'upa-psu');
    tseries = tseries.';
    
    %% build csdm
    
    fprintf('averaging %s...\n',ifilestructs(ixfile).name);
    zphones = (1:24)*.5; zphones = zphones(:)-mean(zphones);
    fs = 102400;
    fftsize = fs/10;
    [csdm1,fr,nsamp1] = Make_CSDMatrix(tseries,fs,fftsize,400,2000,'hanning',.5);
    
    % accumulate between files
    if ~exist('csdm','var')
        csdm=csdm1;
        nsamp=nsamp1;
    else
        csdm=csdm+csdm1;
        nsamp=nsamp+nsamp1;
    end
    
    %% norm csdm
    
    csdm1(:) = Norm_CSDM_Power(csdm,fs,nsamp);
    
    types = { 'conventional', 'mvdr', 'mvdr-symmetric' };
    
    fprintf('beamforming...\n');
    a={};
    for ix=1:length(types)
        type = types{ix};
        % beamform
        [a{ix}.rl,a{ix}.rl_theta,a{ix}.bf,a{ix}.bf_theta] = ...
            Make_RL(csdm1,fr,zphones,100,type,'taylorwin',[]);
        a{ix}.rl = -10*log10(a{ix}.rl);
        a{ix}.bf = 10*log10(a{ix}.bf);
    end
    
    bfmax=-Inf;
    for ix=1:length(types), bfmax=max([a{ix}.bf(:);bfmax]); end
    
    fprintf('plotting...\n');
    clf();
    for ix=1:length(types)
        rl=a{ix}.rl;
        rl_theta=a{ix}.rl_theta;
        bf=a{ix}.bf;
        bf_theta=a{ix}.bf_theta;
        
        %% make plot
        
        %image beamformer output
        subplot(321+2*(ix-1));
        imagesc(fr,bf_theta,bf);caxis([-30 0]+55);
        set(gca,'YDir','normal','YTick',-90:30:90);
        title(sprintf('Beam-former Output (dB re 1 upa)\n%s, %d files',type,ixfile));
        xlabel('Frequency (Hz)');
        ylabel('Beam angle (deg)');
        colorbar;
        
        %image RL
        subplot(322+2*(ix-1));
        imagesc(rl_theta,fr,flipud(rl));caxis([0,10])
        set(gca,'YDir','normal','XTick',0:30:90);
        title(sprintf('Reflection Loss (dB)\n%s, %d files',type,ixfile));
        xlabel('Grazing angle (deg)');
        ylabel('Frequency (Hz)');
        colorbar;
        
    end
    
    drawnow();
    
end

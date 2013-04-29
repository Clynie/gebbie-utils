function [] = Plot_Beamforming(bf,fr,subtitle,font_size,db_range)

if ~exist('subtitle','var'), subtitle = []; end
if ~exist('font_size','var') || isempty(font_size), font_size = 14; end
if ~exist('db_range','var'), db_range = []; end

xx=sort(fr)/1000;
yy=[-90 90];
zz=10*log10(bf);

imagesc(xx,yy,zz);

colorbar();
if ~isempty(db_range), caxis(db_range); end
ylabel('look angle (deg)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
xlabel('frequency (kHz)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
title_str = 'Beam Power (dB re 1 uPa)';
if ~isempty(subtitle), title_str = sprintf('%s\n%s',title_str,subtitle); end
title(title_str,'Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
set(gca,'FontSize',font_size,'FontName','Times New Roman');
set(gca,'YDir','normal');
set(gca,'YTick',-90:15:90);

function [] = Plot_ReflectionLoss(bf,fr,subtitle,font_size,db_range)

if ~exist('subtitle','var'), subtitle = []; end
if ~exist('font_size','var') || isempty(font_size), font_size = 14; end
if ~exist('db_range','var'), db_range = [0 15]; end

nbeams = size(bf,1);

N = floor(nbeams/2);
rl = (bf(nbeams-N+1:end,:)+eps)./(bf(N:-1:1,:)+eps);
rl = rot90(rl);

xx = linspace(0,90,floor(nbeams/2));
yy = sort(fr,'descend')/1000;
zz = 10*log10(rl);

imagesc(xx,yy,zz);

colorbar();
if ~isempty(db_range), caxis(db_range); end
xlabel('grazing angle (deg)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
ylabel('frequency (kHz)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
title_str = 'Bottom Loss (dB)';
if ~isempty(subtitle), title_str = sprintf('%s\n%s',title_str,subtitle); end
title(title_str,'Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
set(gca,'FontSize',font_size,'FontName','Times New Roman');
set(gca,'YDir','normal');
set(gca,'XTick',0:15:90);

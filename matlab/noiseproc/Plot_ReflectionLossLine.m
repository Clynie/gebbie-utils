function [] = Plot_ReflectionLossLine(bf,fr,frplot,subtitle,font_size,db_range)

if ~exist('bf','var'), error('must specify bf'); end
if ~exist('fr','var'), error('must specify fr'); end
if ~exist('frplot','var'), error('must specify frplot'); end
if ~exist('subtitle','var'), subtitle = []; end
if ~exist('font_size','var') || isempty(font_size), font_size = 14; end
if ~exist('db_range','var'), db_range = [-3 12]; end

nbeams = size(bf,1);
nrlbeams = floor(nbeams/2);
nlines = length(frplot);
rl = nan(nrlbeams,nlines);

for ix = 1:nlines
    [~,fr_idx] = min(abs(fr-frplot(ix)));
    rl(:,ix) = (bf(nbeams-nrlbeams+1:end,fr_idx)+eps)./(bf(nrlbeams:-1:1,fr_idx)+eps);
    frplot(ix) = fr(fr_idx);
end

xx = repmat(linspace(0,90,floor(nbeams/2)).',1,nlines);
yy = 10*log10(rl);

plot(xx,yy);

xlabel('grazing angle (deg)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
ylabel('bottom loss (dB)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
title_str = 'Bottom Loss (dB)';
if ~isempty(subtitle), title_str = sprintf('%s\n%s',title_str,subtitle); end
title(title_str,'Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
set(gca,'FontSize',font_size,'FontName','Times New Roman');
set(gca,'YDir','normal');
set(gca,'XTick',0:15:90);
if ~isempty(db_range), ylim(db_range); end

leg = {};
for ix = 1:nlines
    leg{end+1} = sprintf('%d Hz',round(frplot(ix)));
end
hLeg = legend(leg);
set(hLeg,'Location','SouthEast')

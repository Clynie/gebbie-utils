function [] = Plot_Pfath(pfath,time,subtitle,font_size,envelope,yrange,color)

if ~exist('subtitle','var'), subtitle = []; end
if ~exist('font_size','var') || isempty(font_size), font_size = 14; end
if ~exist('envelope','var') || isempty(envelope), envelope = 1; end
if ~exist('yrange','var'), yrange = []; end
if ~exist('color','var'), color = []; end

Nseries = size(pfath,1);
for n=1:Nseries
    pfath(n,:) = pfath(n,:)-mean(pfath(n,:)); % make zero-mean
    if envelope, pfath(n,:) = abs(hilbert(pfath(n,:))); end % envelope, if asked for
end

xx = time.'*1500/2; % 1500=c, 2=two-way travel time
yy = pfath.';
min_idx = min([length(xx) length(yy)]);
hLines = plot(xx(1:min_idx),yy(1:min_idx));

xlabel('depth (m)','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
ylabel('i.r. envelope','Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
title_str = 'Passive Fathometer';
if ~isempty(subtitle), title_str = sprintf('%s\n%s',title_str,subtitle); end
title(title_str,'Interpreter','none','FontSize',font_size,'FontName','Times New Roman');
set(gca,'FontSize',font_size,'FontName','Times New Roman');
if ~isempty(yrange), ylim(yrange); end
xlim([min(xx(:)) max(xx(:))]);
if ~isempty(color)
    if size(color,1)==1 && length(hLines)>1
        color = repmat(color,length(hLines),1);
    end
    for ix=1:length(hLines)
        set(hLines(ix),'Color',color(ix,:));
    end
end

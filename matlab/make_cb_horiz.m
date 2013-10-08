function [hCb] = make_cb_horiz(label, ticks)

if nargin < 1
    label = [];
end

CC = colormap();            % N-by-3
CC = permute(CC, [3 1 2]);  % 1-by-N-by-3

clims = get(gca, 'CLim');

if nargin < 2
    ticks = linspace(clims(1), clims(end), 11);
end

hCb = axes();

clims(1) = min(clims(1), min(ticks));
clims(2) = max(clims(2), max(ticks));

xx = linspace(clims(1), clims(2), size(CC, 2))';
yy = zeros(size(xx));
image(xx, yy, CC);
set(hCb, 'TickDir', 'out');
set(hCb, 'XTick', ticks);
set(hCb, 'YTick', []);
set(hCb, 'box', 'on');
if ~isempty(label)
    xlabel(hCb, label);
end

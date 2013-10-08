function [hCb] = make_cb_vert(label, ticks)

if nargin < 1
    label = [];
end

CC = colormap();            % N-by-3
CC = permute(CC, [1 3 2]);  % N-by-1-by-3

clims = get(gca, 'CLim');

if nargin < 2
    ticks = linspace(clims(1), clims(end), 11);
end

hCb = axes();

clims(1) = min(clims(1), min(ticks));
clims(2) = max(clims(2), max(ticks));

yy = linspace(clims(1), clims(2), size(CC, 1))';
xx = zeros(size(yy));
image(xx, yy, CC);
set(hCb, 'TickDir', 'out');
set(hCb, 'YTick', ticks);
set(hCb, 'XTick', []);
set(hCb, 'box', 'on');
if ~isempty(label)
    xlabel(hCb, label);
end

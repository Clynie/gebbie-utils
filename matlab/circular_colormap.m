function [cmap] = circular_colormap(N, colors)
% Circular colormap.  Useful for coloring angles.
%
% INPUT N : number of steps
%
% INPUT colors : an M-by-3 matrix of colors.  If the last color matches the
% first, the colormap will be circular.
%
% OUTPUT cmap : a colormap that can be passed to the 'colormap' command.
%
% Author: John Gebbie
% Creation Date: 2013 Aug 18

if nargin < 1
    N = 32;
end
if nargin < 2
    colors = [...
        0 0 1 ;...
        0 1 0 ;...
        1 0 0 ;...
        1 1 0 ;...
        0 0 1 ;...
        ];
end

cmap = linear_colormap(N, colors);

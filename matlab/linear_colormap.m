function [cmap] = linear_colormap(N, colors)
% Arbitrary colormap that linearly interpolates between colors spaced
% evenly on the color axis.  Defaults to a brightened jet colormap.
%
% INPUT N : number of steps
%
% INPUT colors : an M-by-3 matrix of colors.  If the last color matches the
% first, the colormap will be circular.
%
% OUTPUT cmap : an N-by-3 matrix of colors. This can be passed directly to
% the colormap command.
%
% Author: John Gebbie
% Creation Date: 2013 Aug 18

if nargin < 1
    N = 32;
end
if nargin < 2
    colors = [...
        0 0 1 ;...
        0 1 1 ;...
        1 1 0 ;...
        1 0 0 ;...
        ];
end

Nc = size(colors,1);        % number of colors

cmap = zeros(N,3);          % initialize output

x = (0:N-1)'/N;             % x-axis from 0 to 1
cx = (0:Nc-1)'/(Nc-1);      % x-value of pure colors
S = 1/(Nc-1);               % spacing between colors on x-axis

for n = 1:Nc
    d = abs(x-cx(n));       % distance from q
    y = 1-d/S;              % strength of color inv. prop. to dist
    y(y<0) = 0;             % define only for region around curr. color
    cmap = cmap + y*colors(n,:);
end

cmap(cmap>1) = 1;           % fix rounding

function [hAxes] = format_axes_v4(cfg,myRows,myCols,fmtWnd,hAxes)
% FORMAT_AXES_V4 Create axes at an exact place on the current figure.
%   INPUTS
%       cfg         => configuration struct (see below for spec)
%       myRows      => vector of row numbers to span
%       myCols      => vector of column numbers to span
%       fmtWnd      => whether to set size of window (boolean). can
%                       specify true first time this is called and false
%                       after to speed up plotting.
%       hAxes       => existing axes to modify
%
%   OUTPUTS
%       hAxes       => handle to new axes
%
%   CFG struct
%       units       => string indicating units of window dimensions
%                       Usualy either 'Inches' or 'Pixels'.
%       width       => width  of entire figure (specified in 'units')
%       height      => height of entire figure (specified in 'units')
%       left        => size of left   gutter (in normalized window units)
%       right       => size of right  gutter
%       top         => size of top    gutter
%       bottom      => size of bottom gutter
%       row         => size of between-row gutter
%       col         => size of between-column gutter
%       nCols       => total number of columns
%       nRows       => total number of rows
%
%   EXAMPLE
%       cfg.units = 'Inches';
%       cfg.width = 3.5;
%       cfg.height = 4;
%       cfg.left = 0.10;
%       cfg.right = 0.10;
%       cfg.top = 0.10;
%       cfg.bottom = 0.10;
%       cfg.row = 0.10;
%       cfg.col = 0.10;
%       cfg.nRows = 1;
%       cfg.nCols = 1;
%       format_axes_v3(cfg,1,1);
%
% Author: John Gebbie
% Modify Date: March 26, 2012

if ~exist('fmtWnd','var') || isempty(fmtWnd)
    fmtWnd = true;
end

% move figure onto the main screen if not already docked
if fmtWnd && ~strcmp(get(gcf,'WindowStyle'),'docked')
    set(gcf,'WindowStyle','normal');
    set(gcf,'Units','Normalized');
    set(gcf,'Position',[.2 .2 .8 .8]);
end

% set the correct units
if ~isfield(cfg,'units') || isempty(cfg.units), cfg.units = 'Inches'; end
set(gcf,'Units',cfg.units);
if strcmpi(cfg.units,'pixels')
    set(gcf,'PaperUnits','points');
else
    set(gcf,'PaperUnits',cfg.units);
end

% set the correct dimensions
if fmtWnd
    set(gcf,'PaperSize',[cfg.width cfg.height]);
    set(gcf,'PaperPosition',[0 0 cfg.width cfg.height]);
    if ~strcmp(get(gcf,'WindowStyle'),'docked')
        position = get(gcf,'Position');
        set(gcf,'Position',[position(1:2) cfg.width cfg.height]);
    end
end

if ~exist( 'hAxes', 'var' ) || isempty( hAxes )
    hAxes = axes();
end

set(hAxes,'Units','normalized');

colSpan         = max(myCols) - min(myCols) + 1;
rowSpan         = max(myRows) - min(myRows) + 1;

colWidth        = ( 1 - cfg.left   - cfg.right - cfg.col*(cfg.nCols-1) ) / cfg.nCols;
rowHeight       = ( 1 - cfg.bottom - cfg.top   - cfg.row*(cfg.nRows-1) ) / cfg.nRows;

axesWidth       = colWidth  * colSpan + cfg.col*(colSpan-1);
axesHeight      = rowHeight * rowSpan + cfg.row*(rowSpan-1);

leftPos         = cfg.left   + (min(myCols)-1)        *(colWidth +cfg.col);
bottomPos       = cfg.bottom + (cfg.nRows-max(myRows))*(rowHeight+cfg.row);

set(hAxes,'Position',[leftPos bottomPos axesWidth axesHeight]);

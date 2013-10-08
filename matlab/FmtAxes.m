classdef FmtAxes < handle
    % Utility for creating multiple axes on a single figure.  This provides
    % greater control over the position of each axis than subplot.  Like
    % subplot, this utility lays out axes in a grid.  Each axes can span
    % multiple rows and/or columns.  Figure dimensions are specified in any
    % units.  Gutters are defined in normalized figure units.
    %
    %   EXAMPLE
    %       o = FmtAxes();
    %       o.units = 'Inches';
    %       o.width = 3.5;
    %       o.height = 4;
    %       o.left = 0.10;
    %       o.right = 0.10;
    %       o.top = 0.10;
    %       o.bottom = 0.10;
    %       o.row = 0.10;
    %       o.col = 0.10;
    %       o.nRows = 1;
    %       o.nCols = 1;
    %       o.fmt();
    %
    % Author: John Gebbie
    % Modify Date: 2013 May 7
    
    properties
        
        % string indicating units of window dimensions. Usualy either
        % 'Inches' or 'Pixels'
        units         = 'Inches'
        
        % width  of entire figure (specified in 'units')
        width         = 3.5
        
        % height of entire figure (specified in 'units')
        height        = 4
        
        % size of left   gutter (in normalized window units)
        left          = 0.10
        
        % size of right  gutter
        right         = 0.10
        
        % size of top    gutter
        top           = 0.10
        
        % size of bottom gutter
        bottom        = 0.10
        
        % size of between-row gutter
        row           = 0.10
        
        % size of between-column gutter
        col           = 0.10
        
        % total number of rows
        nRows         = 1
        
        % total number of columns
        nCols         = 1
        
        % whether to set size of window (boolean). specify true first time
        % axes are created and will automatically switch to false after to
        % speed up plotting.
        fmtWnd        = true
        
        % graphics object to place new axes relative to. defaults to the
        % figure itself. useful for placing relative to another axes in the
        % figure.
        relative_objs  = []
    end
    
    methods
        
        function [hAxes] = fmt(o, myRows, myCols, hAxes)
            % Setup axes at exact place on the current figure.
            %
            % INPUT myRows : vector of row numbers to span
            %
            % INPUT myCols : vector of column numbers to span
            %
            % INPUT hAxes : existing axes to modify (optional, will create
            % new axes if unspecified)
            %
            % OUTPUT hAxes : handle to new axes
            
            if nargin < 2
                % if no arguments are passed in, create the first set of axes
                myRows = 1;
                myCols = 1;
                
            elseif nargin == 2 && ishandle(myRows)
                % handle case if just a axes handle is passe din
                hAxes = myRows;
                myRows = 1;
                myCols = 1;
                
            elseif nargin < 4
                % if axes not specified, define the variable as empty
                hAxes = [];
            end
            
            % move figure onto the main screen if not already docked
            if o.fmtWnd && ~strcmp(get(gcf, 'WindowStyle'), 'docked')
                set(gcf, 'WindowStyle', 'normal');
                set(gcf, 'Units', 'Normalized');
                set(gcf, 'Position', [.2 .2 .8 .8]);
            end
            
            % set the correct units
            set(gcf, 'Units', o.units);
            if strcmpi(o.units, 'pixels')
                set(gcf, 'PaperUnits', 'points');
            else
                set(gcf, 'PaperUnits', o.units);
            end
            
            % set the correct dimensions
            if o.fmtWnd
                set(gcf, 'PaperSize', [o.width o.height]);
                set(gcf, 'PaperPosition', [0 0 o.width o.height]);
                if ~strcmp(get(gcf, 'WindowStyle'), 'docked')
                    position = get(gcf, 'Position');
                    set(gcf, 'Position', [position(1:2) o.width o.height]);
                end
                o.fmtWnd = false;
            end
            
            % create axes, if none passed in
            if isempty(hAxes)
                hAxes = axes();
            end
            
            % ensure axis is in normalized units
            set(hAxes, 'Units', 'normalized');
            
            % compute position of axes on figure canvas
            colSpan         = max(myCols) - min(myCols) + 1;
            rowSpan         = max(myRows) - min(myRows) + 1;
            colWidth        = ( 1 - o.left   - o.right - o.col*(o.nCols-1) ) / o.nCols;
            rowHeight       = ( 1 - o.bottom - o.top   - o.row*(o.nRows-1) ) / o.nRows;
            axesWidth       = colWidth  * colSpan + o.col*(colSpan-1);
            axesHeight      = rowHeight * rowSpan + o.row*(rowSpan-1);
            leftPos         = o.left   + (min(myCols)-1)      *(colWidth +o.col);
            bottomPos       = o.bottom + (o.nRows-max(myRows))*(rowHeight+o.row);
            
            if ~isempty(o.relative_objs)
                
                opos = nan(length(o.relative_objs),4);
                for n = 1:length(o.relative_objs)
                    ro = o.relative_objs(n);
                    ru = get(ro, 'Units');
                    set(ro, 'Units', 'Normalized');
                    opos1 = get(ro, 'Position');
                    set(ro, 'Units', ru);
                    assert(length(opos1) == 4);
                    opos(n,:) = opos1;
                end
                
                oxLL = opos(:,1);
                oyLL = opos(:,2);
                oxUR = oxLL + opos(:,3);
                oyUR = oyLL + opos(:,4);
                ox = [min(oxLL) max(oxUR)];
                oy = [min(oyLL) max(oyUR)];
                
                x = [0; axesWidth] + leftPos;
                y = [0; axesHeight] + bottomPos;
                x = x .* diff(ox) + ox(1);
                y = y .* diff(oy) + oy(1);
                
                leftPos = x(1);
                bottomPos = y(1);
                axesWidth = diff(x);
                axesHeight = diff(y);
                
            end
            
            % put the axes in the calculated place.
            set(hAxes, 'Position', [leftPos bottomPos axesWidth axesHeight]);
        end
        
        
        function [fmt] = duplicate(o)
            % Create duplicate of this class, but using a different handle.
            % Use this to create a copy such that the original object is
            % kept separate and unmodified.
            %
            % OUTPUT fmt : a new FmtAxes object.
            
            fmt = FmtAxes();
            fns = properties(o);
            for ix = 1:length(fns)
                fmt.(fns{ix}) = o.(fns{ix});
            end
        end
        
        
        function [] = set_jasa_width(o, num_cols)
            % Sets the 'width' property for a figure to be included in a
            % JASA article.
            %
            % INPUT num_cols : number of columns to span (1 or 2)
            
            if nargin < 2
                num_cols = 1;
            end
            
            JASA_PAGE_WIDTH = 8+1/2;
            JASA_MARGIN_SIDE = 3/4;
            JASA_COL_WIDTH = 4+1/8 - JASA_MARGIN_SIDE;
            
            if num_cols == 1
                o.width = JASA_COL_WIDTH;
            else
                o.width = JASA_PAGE_WIDTH - 2*JASA_MARGIN_SIDE;
            end
        end
        
    end
    
    methods (Static)
        
        function [] = print(filepath, dpi)
            % Create an image file for the current figure.
            %
            % INPUT filepath : path to file to create. the extension
            % determines the file type. (e.g. png)
            %
            % INPUT dpi : (optional) resolution, defaults to 150.
            
            if nargin < 2
                dpi = [];
            end
            
            % start building argument list, start with handle
            args = { gcf };
            
            % use the file extenstion to determine the rendering engine
            [~, ~, type] = fileparts(filepath);
            type = lower(type(2:end)); % hack off period, and make LC
            switch type
                case ''
                    type = 'png';
                case 'jpg'
                    type = 'jpeg';
            end
            args{end+1} = sprintf('-d%s', type);
            
            % append the resolution, if specified
            if ~isempty(dpi)
                args{end+1} = sprintf('-r%d', dpi);
            end
            
            % append the filepath last
            args{end+1} = filepath;
            
            % call print
            print(args{:});
            
        end
        
        function [] = set_fonts_presentation()
            FmtAxes.set_fonts('Arial', 16, 'normal');
        end
        
        function [] = set_fonts_document()
            FmtAxes.set_fonts('Times New Roman', 12, 'normal');
        end
        
        function [] = set_fonts(fig, font, size, weight, units)
            
            if nargin < 1 || isempty(fig)
                fig = gcf;
            end
            if nargin < 2 || isempty(font)
                font = 'Times New Roman';
            end
            if nargin < 3 || isempty(size)
                size = 9;
            end
            if nargin < 4 || isempty(weight)
                weight = 'normal';
                % options: [ light | {normal} | demi | bold ]
            end
            if nargin < 5 || isempty(units)
                units = 'points';
                % options: [ inches | centimeters | normalized | {points} | pixels ]
            end
            
            handles = findall(fig);
            for n = 1:length(handles)
                
                try
                    set(handles(n), 'FontName', font);
                catch e %#ok
                end
                
                try
                    set(handles(n), 'FontUnits', units);
                catch e %#ok
                end
                
                try
                    set(handles(n), 'FontSize', size);
                catch e %#ok
                end
                
                try
                    set(handles(n), 'FontWeight', weight);
                catch e %#ok
                end
                
            end
        end
        
        
    end
end


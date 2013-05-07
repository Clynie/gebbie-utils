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
    end
    
    methods
        
        function [hAxes] = fmt(o, myRows, myCols, hAxes)
            % FMT Setup axes at exact place on the current figure.
            %
            %   INPUTS
            %       o           => this object
            %       myRows      => vector of row numbers to span
            %       myCols      => vector of column numbers to span
            %       hAxes       => existing axes to modify (optional, will
            %                       create new axes if unspecified)
            %
            %   OUTPUTS
            %       hAxes       => handle to new axes
            
            % if no arguments are passed in, create the first set of axes
            if nargin < 2
                myRows = 1;
                myCols = 1;
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
            
            % create axes, if unspecified
            if nargin < 4 || isempty(hAxes)
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
            
            % put the axes in the calculated place.
            set(hAxes, 'Position', [leftPos bottomPos axesWidth axesHeight]);
        end
    end
    
    methods (Static)
        
        function [] = print(filepath, dpi)
            % Create an image file for the current figure.
            %
            %   INPUTS
            %       o           => this object
            %       filepath    => path to file to create. the extension
            %                       determines the file type. (e.g. png)
            %       dpi         => (optional) resolution, defaults to 150.
            
            if nargin < 2
                dpi = [];
            end
            
            % start building argument list, start with handle
            args = { gcf };
            
            % use the file extenstion to determine the rendering engine
            [~, ~, type] = fileparts(filepath);
            type = lower(type(2:end)); % hack of period, and make LC
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
        
        function [] = set_fonts(font_name, font_size, font_weight)
            h = findall(gcf);
            for n = 1:length(h)
                try
                    set(h(n), 'FontName', font_name);
                catch e %#ok
                end
                try
                    set(h(n), 'FontSize', font_size);
                catch e %#ok
                end
                try
                    set(h(n), 'FontWeight', font_weight);
                catch e %#ok
                end
            end
        end
        
    end
end


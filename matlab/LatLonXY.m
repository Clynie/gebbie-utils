% This class is a facility for converting lat/lon coordinates to and from
% planar XY coordinates.
%
% EXAMPLE:
%
%
% clear;
% clc;
%
% latlon_wash_mon     = [38.889460  -77.035243]; % washington monument
% latlon_whitehouse   = [38.897446  -77.036556]; % whitehouse lawn
% latlon_lincoln      = [38.889271  -77.049827]; % lincoln memorial
% latlon_capitol      = [38.889766  -77.009002]; % capitol building
% latlon_pentagon     = [38.871024  -77.055832]; % pentagon building
% laglon_psu          = [45.509310 -122.681292]; % portland state FAB
%
% obj = LatLonXY();
% obj.set_latlon_reference(latlon_wash_mon);
% coord_latlon = [ ...
%     latlon_whitehouse ; ...
%     latlon_lincoln ; ...
%     latlon_capitol ; ...
%     latlon_pentagon ; ...
%     laglon_psu ];
% coord_xy = obj.convert_latlon_to_xy(coord_latlon);
%
% str_units = '(x,y) (m)';
% fprintf('%20s %10s: (%15f, %15f)\n', 'whitehouse', str_units, coord_xy(1,:));
% fprintf('%20s %10s: (%15f, %15f)\n', 'lincoln memorial', str_units, coord_xy(2,:));
% fprintf('%20s %10s: (%15f, %15f)\n', 'capitol rotunda', str_units, coord_xy(3,:));
% fprintf('%20s %10s: (%15f, %15f)\n', 'pentagon', str_units, coord_xy(4,:));
% fprintf('%20s %10s: (%15f, %15f)\n', 'psu FAB', str_units, coord_xy(5,:));
%
% coord_latlon2 = obj.convert_xy_to_latlon(coord_xy);
% residual = coord_latlon(:) - coord_latlon2(:);
% fprintf('%20s %10s: %e\n', 'residual', 'deg', mean(abs(residual)));
%
%
% Author: John Gebbie
% Creation Date: Feb 7, 2012

classdef LatLonXY < handle
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
        % center of lat/lon grid
        latlon_reference
        
        % raduis of earth
        rEarth = 6378.1*1000; % m
        
        % transform from XYZ on earth sphere to top of sphere making
        % 'flattening' easy
        rotate_earth_xform_mtx
        
    end % properties
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        
        function [] = set_latlon_reference(obj, latlon_reference)
            % [] = set_latlon_reference(obj, latlon_reference)
            %
            %   Set the lat/lon coordinate which will be at the origin of
            %   the XY grid.
            %
            % INPUTS
            %           obj                 : this object
            %           latlon_reference    : a 2-element vector [lat lon]
            %                                   specified in fractional
            %                                   degrees.
            
            
            %//////////////////////////////////////////////////////////////
            % store the center of grid
            %//////////////////////////////////////////////////////////////
            
            obj.latlon_reference = latlon_reference;
            
            %//////////////////////////////////////////////////////////////
            % rotate the reference point, which is on the earth's surface,
            % about the z-axis so that they land directly over the negative
            % y-axis. use negative y-axis instead of positive y-axis so
            % that after the next rotation we still preserve the convention
            % that north points towards the positive y-axis. see:
            % http://en.wikipedia.org/wiki/Rotation_matrix
            %//////////////////////////////////////////////////////////////
            
            ThetaZ = -obj.latlon_reference(2)*pi/180-pi/2;
            Rz = [cos(ThetaZ) -sin(ThetaZ) 0 ; sin(ThetaZ) cos(ThetaZ) 0 ; 0 0 1];
            
            %//////////////////////////////////////////////////////////////
            % rotate again about the x-axis to put the reference point on
            % the very top of the globe, colinear with the z-axis. after
            % this rotation, we then just discard the z coordinate of each
            % transformed point (set it to zero) leaving a flatted 2-D grid
            % of XY points. note that points with more northern latitudes
            % maintain an orientation that points along the positive
            % y-axis.
            %//////////////////////////////////////////////////////////////
            
            ThetaX = -(pi/2-obj.latlon_reference(1)*pi/180);
            Rx = [1 0 0 ; 0 cos(ThetaX) -sin(ThetaX) ; 0 sin(ThetaX) cos(ThetaX)];
            
            %//////////////////////////////////////////////////////////////
            % combine both rotation matrices into a single transform
            % matrix.
            %//////////////////////////////////////////////////////////////
            
            obj.rotate_earth_xform_mtx = Rx * Rz;
            
        end % function
        
        
        function [xy_m] = convert_latlon_to_xy(obj, latlon)
            % [xy_m] = convert_latlon_to_xy(obj, latlon)
            %
            %   Convert lat/lon to XY in meters by rotation and flattening.
            %   Note that set_latlon_reference() must be called prior to
            %   this method.
            %
            % INPUTS
            %           obj         : this object
            %           latlon      : latitude values in first column and
            %                           longitude values in second column.
            %
            % OUTPUT
            %           xy_m        : x and y coordinates of each input
            %                           value specified in meters. the
            %                           first column is for x and the
            %                           second is for y.
            
            %//////////////////////////////////////////////////////////////
            % compute 3-D coordinate of each lat/lon coordinate
            %//////////////////////////////////////////////////////////////
            
            lat = latlon(:, 1);
            lon = latlon(:, 2);
            
            xyz_m = nan(1, numel(lat));
            ThetaLons = lon(:).'*pi/180;
            ThetaLats = lat(:).'*pi/180;
            xyz_m(1, :) = obj.rEarth .* cos(ThetaLons) .* cos(ThetaLats);
            xyz_m(2, :) = obj.rEarth .* sin(ThetaLons) .* cos(ThetaLats);
            xyz_m(3, :) = obj.rEarth .* sin(ThetaLats);
            
            %//////////////////////////////////////////////////////////////
            % rotate points to be on the top of the sphere
            %//////////////////////////////////////////////////////////////
            
            xyz_m = obj.rotate_earth_xform_mtx * xyz_m;
            
            %//////////////////////////////////////////////////////////////
            % configure output variable
            %//////////////////////////////////////////////////////////////
            
            xy_m = xyz_m(1:2, :).';
            
        end % function
        
        
        function [latlon] = convert_xy_to_latlon(obj, xy)
            % [latlon] = convert_xy_to_latlon(obj, xy)
            %
            %   Convert x/y coordinates to lat/lon coordinates.  Note that
            %   set_latlon_reference() must be called prior to this method.
            %
            % INPUTS
            %           obj         : this object
            %           xy          : xy coordinates where first column
            %                           contains the x components and the
            %                           second column contains the y
            %                           components.
            %
            % OUTPUTS
            %           latlon      : latitude and longitude coordinates
            %                           for each input XY coordinate. the
            %                           first column is for latitudes and
            %                           the second column is for
            %                           longitudes.
            
            %//////////////////////////////////////////////////////////////
            % re-construct Z coordinate from x and y
            %//////////////////////////////////////////////////////////////
            
            x = xy(:, 1);
            y = xy(:, 2);
            
            z = sqrt(obj.rEarth.^2 - x.^2 - y.^2);
            
            %//////////////////////////////////////////////////////////////
            % rotate points to be back at their original place on the
            % sphere
            %//////////////////////////////////////////////////////////////
            
            xyz_m = [x(:).' ; y(:).' ; z(:).'];
            xyz_m = obj.rotate_earth_xform_mtx \ xyz_m;
            
            %//////////////////////////////////////////////////////////////
            % use trig to figure out lat/lon angles on the sphere
            %//////////////////////////////////////////////////////////////
            
            lat = asin(xyz_m(3, :)./obj.rEarth)*180/pi;
            lon = atan2(xyz_m(2, :), xyz_m(1, :))*180/pi;
            
            %//////////////////////////////////////////////////////////////
            % configure output variable
            %//////////////////////////////////////////////////////////////
            
            latlon = [lat(:) lon(:)];
            
        end % function
        
    end % methods
    
end % classdef

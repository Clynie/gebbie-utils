% Class for doing linear transformations to 3-D coordinates.
%
% EXAMPLE:
%
%     % initialize workspace
%     clear;
%     clc;
%
%
%     % Create some sample points
%     sample_pts_local = [...
%         0 0 0 ; ...
%         0.5 0 0 ; ...
%         1 0 0 ; ...
%         1.25 0 0 ; ...
%         1.125 .25 0 ; ...
%         1.125 -.25 0 ; ...
%         1 .5 0 ; ...
%         1 -.5 0 ; ...
%         ]';
%
%     % Create an object for converting from coordinate system to another
%     t = TransformationMatrix();
%
%     deg2rad = @(x) x.*pi./180; % function to covert degrees to radians
%
%     % rotate about the z axis (angle is in radians)
%     t.rotate_z(deg2rad(45));
%
%     % rotate about the x axis (angle is in radians)
%     t.rotate_x(deg2rad(60));
%
%     % move to another point in space
%     t.translate(3, -2, 4);
%
%     % apply the transformation
%     sample_pts_global = t.apply(sample_pts_local);
%
%     % clear the plot
%     clf;
%
%     % plot the local points in blue
%     scatter3(sample_pts_local(1,:),sample_pts_local(2,:),sample_pts_local(3,:),'bo');
%     hold on;
%
%     % plot the global points in red
%     scatter3(sample_pts_global(1,:),sample_pts_global(2,:),sample_pts_global(3,:),'ro');
%
%     % setup axes
%     axis tight;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%
%
% Author: John Gebbie
% Creation Date: Nov 11, 2012

classdef TransformationMatrix < handle
    properties
        M = eye(4)          % transformation matrix
        Mi = []             % inverse transformation matrix
    end
    
    methods
        function [] = translate(o,dx,dy,dz)
            % move by [dx,dy,dz]
            A = [ 1 0 0 dx ; 0 1 0 dy ; 0 0 1 dz ; 0 0 0 1 ];
            o.M = A * o.M;
            o.Mi = [];
        end
        
        function [] = reflect_x(o,x)
            o.translate(-x,0,0);
            o.scale_dim(-1,1,1);
            o.translate(x,0,0);
        end
        
        function [] = reflect_y(o,y)
            o.translate(0,-y,0);
            o.scale_dim(1,-1,1);
            o.translate(0,y,0);
        end
        
        function [] = reflect_z(o,z)
            o.translate(0,0,-z);
            o.scale_dim(1,1,-1);
            o.translate(0,0,z);
        end
        
        function [] = rotate_x(o,theta_rad,about_pt)
            % counter-clockwise when the axis of rotation points toward the
            % observer, and the coordinate system is right-handed
            if nargin < 3
                about_pt = [];
            end
            if ~isempty(about_pt)
                o.translate(-about_pt(1),-about_pt(2),-about_pt(3));
            end
            c = cos(theta_rad);
            s = sin(theta_rad);
            A = [ 1 0 0 0 ; 0 c -s 0 ; 0 s c 0 ; 0 0 0 1 ];
            o.M = A * o.M;
            if ~isempty(about_pt)
                o.translate(about_pt(1),about_pt(2),about_pt(3));
            end
            o.Mi = [];
        end
        
        function [] = rotate_y(o,theta_rad,about_pt)
            % counter-clockwise when the axis of rotation points toward the
            % observer, and the coordinate system is right-handed
            if nargin < 3
                about_pt = [];
            end
            if ~isempty(about_pt)
                o.translate(-about_pt(1),-about_pt(2),-about_pt(3));
            end
            c = cos(theta_rad);
            s = sin(theta_rad);
            A = [ c 0 s 0 ; 0 1 0 0 ; -s 0 c 0 ; 0 0 0 1 ];
            o.M = A * o.M;
            if ~isempty(about_pt)
                o.translate(about_pt(1),about_pt(2),about_pt(3));
            end
            o.Mi = [];
        end
        
        function [] = rotate_z(o,theta_rad,about_pt)
            % counter-clockwise when the axis of rotation points toward the
            % observer, and the coordinate system is right-handed
            if nargin < 3
                about_pt = [];
            end
            if ~isempty(about_pt)
                o.translate(-about_pt(1),-about_pt(2),-about_pt(3));
            end
            c = cos(theta_rad);
            s = sin(theta_rad);
            A = [ c -s 0 0 ; s c 0 0 ; 0 0 1 0 ; 0 0 0 1 ];
            o.M = A * o.M;
            if ~isempty(about_pt)
                o.translate(about_pt(1),about_pt(2),about_pt(3));
            end
            o.Mi = [];
        end
        
        function [] = scale_dim(o,sx,sy,sz)
            % scale each dimension separately
            A = [ sx 0 0 0 ; 0 sy 0 0 ; 0 0 sz 0 ; 0 0 0 1 ];
            o.M = A * o.M;
            o.Mi = [];
        end
        
        function [] = scale(o,s)
            % scale each dimension by s
            o.scale_dim(s,s,s);
            o.Mi = [];
        end
        
        function [p] = apply(o,p)
            % applies transformation to points p in which each point is a
            % column vector
            N = size(p,2);
            p = [ p ; ones(1,N) ];
            p = o.M * p;
            p = p(1:3,:) ./ repmat(p(4,:),3,1);
        end
        
        function [p] = apply_inverse(o,p)
            % applies inverse transformation to points p in which each
            % point is a column vector
            o.check_inverse_computed();
            N = size(p,2);
            p = [ p ; ones(1,N) ];
            p = o.Mi * p;
            p = p(1:3,:) ./ repmat(p(4,:),3,1);
        end
        
        function [] = check_inverse_computed(o)
            % check the inverse is computed
            if isempty(o.Mi)
                o.Mi = inv(o.M);
            end
        end
    end
    
    methods (Static)
        function [] = test()
            p = [1 1 1]';
            o = TransformationMatrix();
            o.rotate_x(pi/4,[0 0 10]');
            p = o.apply(p);
            display(p);
        end
    end
end

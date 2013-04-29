% Example code demonstrating how to use the TransformationMatrix class to
% perform linear transformations.
%
% Author: John Gebbie
% Creation Date: 27 Dec 2012


% initialize workspace
clear;
clc;


% Create some sample points
sample_pts_local = [...
    0 0 0 ; ...
    0.5 0 0 ; ...
    1 0 0 ; ...
    1.25 0 0 ; ...
    1.125 .25 0 ; ...
    1.125 -.25 0 ; ...
    1 .5 0 ; ...
    1 -.5 0 ; ...
    ]';

% Create an object for converting from coordinate system to another
t = TransformationMatrix();

deg2rad = @(x) x.*pi./180; % function to covert degrees to radians

% rotate about the z axis (angle is in radians)
t.rotate_z(deg2rad(45));

% rotate about the x axis (angle is in radians)
t.rotate_x(deg2rad(60));

% move to another point in space
t.translate(3, -2, 4);

% apply the transformation
sample_pts_global = t.apply(sample_pts_local);

% clear the plot
clf;

% plot the local points in blue
scatter3(sample_pts_local(1,:),sample_pts_local(2,:),sample_pts_local(3,:),'bo');
hold on;

% plot the global points in red
scatter3(sample_pts_global(1,:),sample_pts_global(2,:),sample_pts_global(3,:),'ro');

% setup axes
axis tight;
xlabel('x');
ylabel('y');
zlabel('z');


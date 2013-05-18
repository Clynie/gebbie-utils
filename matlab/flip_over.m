function [xyz] = flip_over(xyz, boundary_z)
% Reflect points over XY planes at different points along the Z axis.

for ii = 1:size(boundary_z, 2)
    z = xyz(:, 3)-boundary_z(:, ii);
    z = boundary_z(:, ii) - z;
    xyz = [xyz(:, 1:2) z];
end

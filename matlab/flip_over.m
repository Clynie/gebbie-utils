function [zxy] = flip_over(zxy,boundary_z)
for ii = 1:size(boundary_z,2)
    z = zxy(:,1)-boundary_z(:,ii);
    z = boundary_z(:,ii) - z;
    zxy = [z zxy(:,2:3)];
end

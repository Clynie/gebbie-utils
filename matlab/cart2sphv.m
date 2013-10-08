function [ae] = cart2sphv(xyz)

[az,el] = cart2sph(xyz(:,1), xyz(:,2), xyz(:,3));
ae = [az,el];

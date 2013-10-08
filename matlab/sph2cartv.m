function [xyz] = sph2cartv(ae)

[x,y,z] = sph2cart(ae(:,1), ae(:,2), 1);
xyz = [x,y,z];

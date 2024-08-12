function [mass, distance] = fluxGrid(lon,lat,lon_int,lat_int,depth,d_int,rho,SurfRadius,detector)

% -- CALCULATE MASS OF VOXELS --
% - Calculate bounds of voxel top (so the surface of the voxel) -
latbot = lat  - (lat_int/2); lattop = lat  + (lat_int/2);
lonleft = lon - (lon_int/2); lonright = lon + (lon_int/2); 

% Convert to Radians
latbot = (latbot + 90) * pi/180; % pi/180 converts degree to radian
lattop = (lattop + 90) * pi/180;
lonleft = (lonleft + 180) * pi/180;
lonright = (lonright + 180) * pi/180;

r1 = SurfRadius - (depth + (0.5.*d_int)); % top Radius (m)
r2 = SurfRadius - (depth - (0.5.*d_int)); % bottom Radius (m)

a = 1/3;
vol = ((a*r2.^3) - (a*r1.^3)).*(-cos(lattop) + cos(latbot)).*(lonright-lonleft); %m3
mass = vol.*rho; %kg


if istable(detector) == 1
d1 = table2array(detector(1,1));
d2 = table2array(detector(1,2));
d3 = table2array(detector(1,3));

else
d1 = detector(1,1);
d2 = detector(1,2);
d3 = detector(1,3);

end

% Edited by LS 2021 to stop throwing Undefined function 'uminus' for input arguments of type 'table'.
% An UndefinedFunction error was thrown on the workers for 'uminus'.
% This might be because the file containing 'uminus' is not accessible
% on the workers.  Use addAttachedFiles(pool, files) to specify the
% required files to be attached.  For more information see the
% documentation for 'parallel.Pool/addAttachedFiles'.


% -- CALCULATE DISTANCE TO DETECTOR --
cart_det = cart(d1,d2,d3); %detector(3) is already radius

cart_cell = cart(lon,lat,SurfRadius-depth);

% - Calculate distance - 
%   You are just finding the hyponuse but in 3-D
distance = sqrt((cart_det(1)-cart_cell(:,1)).^2 + (cart_det(2)-cart_cell(:,2)).^2 + ...
                (cart_det(3)-cart_cell(:,3)).^2); %m

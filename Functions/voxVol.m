function vol = voxVol(surfRadius,depth,thickness,lat1,lat2,lon1,lon2)
% VOXVOL calculates the volume of a voxel bounded by lat1, lat2, lon1, and
% lon2, which are IN DEGREES.  surfRadius (m) is the radius to the surface
% of the earth (at that cell), depth (m) is the depth to the center of the
% voxel, thickness (m) is the thickness of the voxel, and rho (kg/m3) is the density. 
%
%   volume = voxVol(surfRadius,depth,thickness,lat1,lat2,lon1,lon2)
%
%
%   Code:
%       a = 1/3;
%       vol = ((a*r2.^3) - (a*r1.^3)).*(-cos(lat2) + cos(lat1)).*(lon2-lon1);
% 
%       volume equation is the triple integral of the sphere equation 
%           p^2*sin(phi) d(p)d(phi)d(theta) over the limits of latitude (phi),
%           longitude (theta), and radius (p)
%
%   UNITS:
%   surfRadius = meters
%   depth = meters
%   thickness = meters
%   lat1,2 = degree
%   lon1,2 = degree
%
%
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        June, 2016                      ----- 
%
% 

lat1 = (lat1 + 90) * pi/180;
lat2 = (lat2 + 90) * pi/180;
lon1 = (lon1 + 180) * pi/180;
lon2 = (lon2 + 180) * pi/180;


r1 = surfRadius - (depth + (0.5.*thickness)); % Bottom Radius (m)
r2 = surfRadius - (depth - (0.5.*thickness)); % Top Radius (m)

a = 1/3;
vol = ((a*r2.^3) - (a*r1.^3)).*(-cos(lat2) + cos(lat1)).*(lon2-lon1); %m3

function [distance] = dis(lon, lat, detector, surfRadius, depth)
%
%DIS calculates the absolute straight distance between two points in meters. 
%Two points are defined by longitude, latitude, and radius.  For points at
%depth (i.e. middle crust), the radius should be equal to the radius of the
%Earth + surface elevation - depth to middle of layer. Longitude (lon) in degrees,
%latitude in degrees, detector = 1x3 array with [lon, lat, radius], surfRadius
%should be radius at that point (r Earth + elevation), and depth is depth to
%center of layer (m)(not center of layer with detector). Set depth = 0 if 
%you want distance to top of crustal column. 
%
% [distance] = dis(longitude, latitude, [3x1 det coords], surfRadius, depth
% to layer)
%
%   longitude = degree
%   latitude = degree
%   det coords = degree lon, degree lat, meters radius
%   surfRadius = meters
%   depth = meters (positive value if deep in the Earth)
%
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        August, 2016                      ----- 
%
%See also dist, distance
%

if istable(detector)==1
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


% - Convert to cartesian coordinates (3-dimensions) - 
cart_det = cart(d1,d2,d3); %detector(3) is already radius

cart_cell = cart(lon,lat,surfRadius-depth);

% - Calculate distance - 
%   You are just finding the hyponuse but in 3-D
distance = sqrt((cart_det(1)-cart_cell(:,1)).^2 + (cart_det(2)-cart_cell(:,2)).^2 + ...
                (cart_det(3)-cart_cell(:,3)).^2); %m




function [lon2,lat2,depth2,lon_int2,lat_int2,d_int2] = miniVox(lon,lat,depth,lon_int,lat_int,d_int,DesSize,SurfRadius)

% MINIVOX is a function that calculates the coordinates (latitude,
% longitude, and depth) of cells with size governed by "DesSize". The new
% cells will be smaller than the original cell size.  If "DesSize" is not
% larger than 1.33*depth, than then "depth2" will equal "depth". 
%
%   INPUTS: 
% lat      = (deg) latitude of center of voxel to be gridded
% lon      = (deg) longitude of voxel of voxel to be gridded
% depth    = (m) depth to center of voxel to be gridded
% lon_int  = (deg) current size of voxel for longitude (degree)
% lat_int  = (deg) current size of voxel for latitude (degree)
% d_int    = (m) current thickness of cell (depth interval)
% DesSize  = (m) size of desired voxel (cubic)
% SurfRadius = (m) radius from center of Earth to surface of Earth

%   OUTPUTS:
% lon2     = (deg) new longitude coordinates
% lat2     = (deg) new latitude coordinates
% depth2   = (m) new depths to center of each cell
% lon_int2 = (deg) new size of voxel for longitude
% lat_int2 = (deg) new size of voxel for latitude
% d_int2   = (deg) new size of voxel for depth

%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        March, 2018                      ----- 
%


% - Re-define variable for the desired size of each cell - 
num = DesSize; %(m)

% - Calculate radius to middle of layer - 
radius = SurfRadius - depth; %(m)

%% Grid Longitude


% -- Find size of 1 lon2 degree (m)(Haversine Formula) at center of cell --
% note: Haversine is accurate at small distances, Spherical law of cosines is not
% We are doing this so we can approximately get a voxel the size we want
% (it wont be perfectly what we define)
r = pi/180; 
x = sqrt(cos(lat*r)*cos(lat*r)*sin(((lon_int*r)/2).^2));
deg_size_lon = 2 * radius * atan2(x,1-x); 

    
% -- Calculate new longitude interval --
lon_int2 = lon_int/round(deg_size_lon/num);  % deg_size/num will yield number of splits, then 1/that to get interval          
    

% -- Calculate new longitude coordinates --
    % "lon - lon_int" brings us to side of cell, then "+(lon_int/2)" brings us to center of first cell.

lon2 = lon - (lon_int/2) + (lon_int2/2) : lon_int2 : lon + (lon_int/2) - (lon_int2/2); 

%% Grid Latitude
 % -- Calculate new latitude interval --
        % Calculate radius of Earth and divide by radius of Earth to yield
        % size (m) of 1 degree. Multiply this size by the latitude interval
        % to get size of one cell. This is the latitude distance at the
        % middle of the cell. 
deg_size_lat = 2*pi*radius/360*lat_int; %(m)
    
% -- Calculate new latitude interval --
lat_int2 = lat_int/round(deg_size_lat/num); 

% -- Calculate new latitude coordinates --
lat2 = lat - (lat_int/2) + (lat_int2/2) : lat_int2 : lat + (lat_int/2) - (lat_int2/2);

%% Grid thickness
    
if d_int > 1.3333*num % if not 1.33* larger, then dont split it up
    d_int2 = d_int./ceil(d_int./num); %(m) thickness of each new voxel
    % - deep goes from top to bottom since depth is positive
    depth2 = (depth-(d_int/2))+(d_int2/2) : d_int2 : (depth+(d_int/2))-(d_int2/2); %(m)
else %dont grid depth because its already small
    depth2 = depth; d_int2 = d_int;
end

%% Convert latitude, longitude, and depth to single column variable
% Create grid of values 
[lon2,lat2,depth2] = meshgrid(lon2,lat2,depth2); 

lon2 = lon2(:); %extract into single column
lat2 = lat2(:); 
depth2 = depth2(:);



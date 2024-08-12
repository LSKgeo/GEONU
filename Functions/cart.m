function [output] = cart(deg_lon,deg_lat,radius)
% CART  Converts the spherical coordinates phi (longitude), theta
% (latitude), and radius to cartesian coordinates.  Phi
% and Theta MUST BE IN degrees. Radius is distance from center of Earth in m.
%
%   coordinates = cart(Rlon, Rlat, radius)
%
%       radian = degree * pi/180
%
%       x = radius .* cos(theta) .* sin(phi);
%       y = radius .* sin(theta) .* sin(phi);
%       z = radius .* cos(phi);
%
%     Convert to radians from degrees:
%       phi = (phi - 180).*pi/180;
%       theta = (90 - theta).*pi/180;
%
%  For more information, see 
%<a href="matlab: web('http://astro.uchicago.edu/cosmus/tech/latlong.html')">http://astro.uchicago.edu/cosmus/tech/latlong.html</a>.
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        June, 2016                      ----- 
%
% See also sph2cart 

x = -radius .* cos(deg_lat* pi/180) .* cos(deg_lon* pi/180);
y = radius .* sin(deg_lat* pi/180);
z = radius .* cos(deg_lat* pi/180) .* sin(deg_lon* pi/180);


output = [x y z];
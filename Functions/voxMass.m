function [mass,vol] = voxMass(SurfRadius,depth,thickness,rho,latbot,lattop,lonleft,lonright,iter)
% VOXMASS calculates the mass of a voxel bounded by latbot,lattop,lonleft,
% lonright, which are IN DEGREES and some thickness.  surfRadius is the 
% radius to the surface of the earth (at that cell in m)), depth is the 
% depth to the center of the voxel (m), thickness is the thickness of the 
% voxel (m), and rho is the density (kg/m3). The latitude and longitude
% values are the latitude at the bottom of your cell ('latbot') and top of 
% cells ('lattop'). 
%
%   mass = voxMass(surfRadius,depth,thickness,rho,latbot,lattop,lonleft,lonright)
%
%
%   Code:
%       a = 1/3;
%       vol = ((a*r2.^3) - (a*r1.^3)).*(-cos(lattop) + cos(latbot)).*(lonright-lonleft);
%       mass = vol.*rho;
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        June, 2016                      ----- 
%
% 
latbot = (latbot + 90) * pi/180; % pi/180 converts degree to radian
lattop = (lattop + 90) * pi/180;
lonleft = (lonleft + 180) * pi/180;
lonright = (lonright + 180) * pi/180;


r1 = SurfRadius - (depth + (0.5.*thickness)); % bottom Radius (m)
r2 = SurfRadius - (depth - (0.5.*thickness)); % top Radius (m)

a = 1/3;
%vol = ((a*r2.^3) - (a*r1.^3)).*(-cos(lattop) + cos(latbot)).*(lonright-lonleft); %m3
% writing out the multiplication instead of using power is faster. probably
% due to the matrix multiplication libraries MATLAB is built on
% AI answer: In MATLAB, element-wise operations on arrays can have different performance characteristics depending on the specific operation being performed. When you multiply an array by itself using A .* A, MATLAB performs an element-wise multiplication directly. This operation is straightforward and involves a simple multiplication of corresponding elements in the array.
vol = ((a*r2.*r2.*r2) - (a*r1.*r1.*r1))        .*  repmat((-cos(lattop) + cos(latbot)).*(lonright-lonleft),1,iter)'; %m3
mass = vol.*rho; %kg

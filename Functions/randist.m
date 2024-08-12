function [dist] = randist(mu,error,int,locate)
%RANDIST normal distribution.
%  RANDIST(mu,error,int,locate) creates a normal distribution (gaussian) 
%  from input variables.  mu = mean, error = +- uncertainty, int = size/type
%  of distribution output, and locate = location on PDF.  
%
%   int = 1, dist will be the central value always.
%   int = 0, returned value is single random value within distribution
%           (unless "locate" >1 in length)
%   int > 0, returns distribution comprised of "int" values
%
%   locate = "single value" will return the random number at "locate"
%       location on the distribution. This is used to provide correlated 
%       values betwen creation of different distributions. If locate is a
%       distribution it should be a uniform distribution (i.e. built using
%       the RANDN function, not RAND). 'locate' = 0 will return mu, while 
%       'locate' = 1 will return (5+3) = 8. 
%
%   
%
%   Example usage: 
%   [dist] = randist(3,2,0) returns a single value from the distribution
%       with mu = 3 and uncertainty +- = 2. 
%
%   [dist] = randist(3,2,0,0.9) returns a single value from the higher end
%       of the distribution with mu = 3 and uncertainty +- = 2. 
%
% See also randn, logdist. 
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                    Created June, 2016                  ----- 
%   -----                     Update June, 2018                  ----- 
%



if nargin == 3 
    if int == 1   % returns central value

                dist = mu;


    elseif int > 1% returns distribution of "int" length 

                dist = randn(int,1).*error + mu;


    elseif int == 0 %returns 1 random value from the distribution 

                dist = randn(1,1).*error + mu;


    end
    
elseif nargin == 4 % Find random numbers when specified random distribution
     if int == 1   % returns central value

                dist = mu;


    elseif int > 1% returns distribution of "int" length 

                dist = locate.*error + mu;


    elseif int == 0 %returns 1 random value from the distribution 

                dist = locate.*error + mu;

     end   
else
    warning('You entered %d input(s). You need 3 or 4.',nargin);
end
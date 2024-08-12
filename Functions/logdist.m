function [dist, log_mu, log_sigma] = logdist(mu,plusError,minusError,int,locate)
% LOGDIST Creates a log-normal distribution (non-gaussian) from input
% variables or selects a single value from a distribution. It is also 
% possible to select the location of the random value on the distribution 
% using the 'locate' parameter (where 0 = mean). 
%
% logdist(mu,plusError,minusError, int, locate)
%
%   mu          = mean 
%   plusError   = positive uncertainty
%   minusError  = negative uncertainty
%   int         = size of distribution ("0" if want single value)
%   locate      = location of value on PDF (can be larger than 1)
%   
%   If int = 1, dist will be the central value always. If int = 0, then 
%   output will be a single random value within the distribution, unless 
%   "locate" is larger than 1. In that case the number of random values
%   will be determined by the length of "locate".
%   
%   [dist, log_mu, log_sigma] = logdist(mu,plusError,minusError,int,locate)
%   provides the distribution (dist), the log(mean) of the distribution,
%   and log(sigma) of the distribution.  log(sigma) is a combination of
%   the positive and negative uncertainty and is calculated as: 
%
%   log_plus = log(mu + plusError) - log_mu;
%   log_minus = log_mu - log(mu - minusError);
%   log_sigma = (log_plus+log_minus)/2; 
%   
% See also randist, lognrnd, lognfit, lognpdf, randn.
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                     Created June, 2016                 ----- 
%   -----                     Updated June, 2018                 ----- 




if nargin == 4 

    if int == 1 % CENTRAL Value
              dist = mu;



    elseif int >1% MONTE CARLO 
            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(randn(int,1).*log_sigma(:)+log_mu(:));




    elseif int ==0 %(i.e. you want 1 random value from the distribution)

            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(randn(1,1).*log_sigma(:)+log_mu(:));


    end
    
elseif nargin == 5
    % Find random numbers when specified random distribution
    if int == 1 % CENTRAL Value
              dist = mu;



    elseif int >1% MONTE CARLO 
            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(locate.*log_sigma(:)+log_mu(:)); %single outputs in single precision




    elseif int ==0 %(i.e. you want 1 random value from the distribution)

            log_mu = log(mu);
            log_plus = log(mu + plusError) - log_mu;
            log_minus = log_mu - log(mu - minusError);
            log_sigma = (log_plus+log_minus)/2;

            dist = exp(locate.*log_sigma(:)+log_mu(:));%single outputs in single precision

    else
        warning('You entered %d input(s). You need 4 or 5.',nargin);
    end  
end

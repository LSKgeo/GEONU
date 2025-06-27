function [stat] = stat2(PDF,whatKind)
% STAT characterizes the input distribution 'PDF' which is a single 
% dimension variable (i.e. only 1 column but many rows or opposite).  The 
% type of characterization is determined by 'whatKind' (see below). Note: 
% output is always a 1x3 matrix. If desired statistics only provide two
% numbers (such as [mean,nanstd]) then the 3rd value will be 'NaN'. Output is
% always central value first, followed by the positive uncertainty, then
% the negative uncertainty (only one value if symmetrical uncertainty). 
%
%
%   SYNTAX: 
%   stat(PDF) finds the median and +-68% confidence limits of the 
%       distribution "PDF" 
%
%   stat(PDF, whatKind) returns the median, + sigma, - sigma , of the 
%       distribution "PDF" correlating to whatKind = 1-7 respectively 
%       (see below). 
%
%   Example: 
%       x = randn(100,1)*2 + 10; 
%       stat(x,4) = 9.6 2.2 1.8 (geoMean, + uncertainty, - uncertainty)
%
%   whatKind: 
%       1   |   Median +- 68.24% confidence limit (~1 sigma, but not)
%       2   |   Median +- 95.45% confidence limit (~2 sigma, but not)
%       3   |   Mean with standard deviation (1 sigma)
%       4   |   Geometric mean with +- nanstd of log-distribution 
%       5   |   normal dist. maximum Liklihood fit
%       6   |   log-normal dist. maximum Likelihood fit
%       7   |   gamma dist. maximum liklihood fit
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                        June, 2016                      -----
%   -----                last modified September, 2018           -----
%
%   See also median, mean, nanstd, mle

n = length(PDF(1,:));

% -- If only 1 function argument, i.e find median value --
if nargin == 1
    whatKind = 1;
end


% -- If 2 function argument, i.e find median value, but specify whatKind = [1-7] --
%int = length(PDF(:,1)); %size of PDF
dim = 1;

% find Median, + sigma, - sigma 
if whatKind == 1 % 68.24% of data (~1 sigma)

   % CALL prctile on 16th, 50th, and 84th percentiles
   stat(:,1:3) = prctile(PDF,[50, 84.12, 15.88 ],dim)';% - PDF_median;%;s_PDF(round((int)/2) + round(0.3413*int)) - PDF_median; % plus sigma




elseif whatKind == 2 % 95.45% of data (~2 sigma)
    stat(:,1:3) = prctile(PDF,[50, 97.725, 2.275],dim)';

elseif whatKind == 3 % mean +- nanstd
   stat(:,1) = mean(PDF,dim,'omitnan')'; 
   stat(:,2) = std(PDF,dim,'omitnan')'; 
   stat(:,3) = nan(n,1); 

elseif whatKind == 4 % mean +- from nanstd of naturalLog(pdf)
    if sum(PDF<0) > 0
        PDF(PDF<0)= 0; % CHANGED THIS LS
        error('Data contains negative values! Cannot fit a log-normal distribution.')
    end
    
   stat = zeros(n,3);
   s_PDF = log(PDF);
   x = std(s_PDF,dim,'omitnan'); 
   PDF_mean = exp(mean(s_PDF,dim,'omitnan'));

   stat(:,1) = PDF_mean';
   stat(:,2) = (exp(mean(s_PDF,dim,'omitnan') + x) - PDF_mean)'; % plus sigma
   stat(:,3) = (PDF_mean - exp(mean(s_PDF,dim,'omitnan') - x))'; % minus sigma     

elseif whatKind == 5 % normal maximum likelihood fit
   stat = mle(PDF,'distribution','normal'); 
   stat(:,3) = nan(n,1);
   
elseif whatKind == 6 % log-normal maximum liklihood fit
    if sum(PDF<0) > 0
        error('Data contains negative values! Cannot fit a log-normal distribution.')
    end    
    
   x = mle(PDF,'distribution','lognormal'); 

   stat(1,1) = exp(x(1));
   stat(1,2) = exp(x(1) + x(2)) - exp(x(1)); % plus sigma
   stat(1,3) = exp(x(1)) - exp(x(1) - x(2)); % minus sigma     

elseif whatKind == 7 % gamma maximum likelihood fit
    if sum(PDF<0) > 0
        error('Data contains negative values! Cannot fit a log-normal distribution.')
    end       
  stat = mle(PDF,'distribution','gamma');   
        
else
     % Provide warning and stop program when wrong inputs
 error('Wrong input of "whatKind". See "help stat" for options') 
 end


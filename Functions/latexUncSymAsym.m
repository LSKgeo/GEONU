function [output] = latexUncSymAsym(data,deci,final)
% LATEXUNCSYMASYM converts the data into latex number. The function checks
% if the data has symmetrical uncertainty or asymetric and outputs the
% latex code accordingly. 
%
%   data = 1x3 matrix (3rd column can be 'NaN'
%   deci = either 1 or 2 = number of decimals 
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                    September, 2018                     -----

if round(data(1),2) == 2.7 ||  round(data(1),2) == 10.5 ||  round(data(1),2) == 2.32
    deci = 1;
end


if deci == 1
        if sum(isnan(data))>0 
    output = sprintf(' %0.1f $%s$%0.1f  &',data(1,1),'\pm',data(1,2)); % symetric unc.
        else
    output = sprintf(' %0.1f $ ^{+%0.1f}_{-%0.1f} $ &',data(1,1),data(1,2),data(1,3));  % asymetric unc.
        end
    
elseif deci == 2
        if sum(isnan(data))>0 
    output = sprintf(' %0.2f $%s$%0.2f &',data(1,1),'\pm',data(1,2)); % symetric unc.
        else
    output = sprintf(' %0.2f $ ^{+%0.2f}_{-%0.2f} $ &',data(1,1),data(1,2),data(1,3));  % asymetric unc.
        end
end


if nargin == 3 %dont have '&' at end
    if deci == 1
        if sum(isnan(data))>0 
    output = sprintf(' %0.1f $%s$%0.1f',data(1,1),'\pm',data(1,2)); % symetric unc.
        else
    output = sprintf(' %0.1f $ ^{+%0.1f}_{-%0.1f} $',data(1,1),data(1,2),data(1,3));  % asymetric unc.
        end
    
elseif deci == 2
        if sum(isnan(data))>0 
    output = sprintf(' %0.2f $%s$%0.2f',data(1,1),'\pm',data(1,2)); % symetric unc.
        else
    output = sprintf(' %0.2f $ ^{+%0.2f}_{-%0.2f} $',data(1,1),data(1,2),data(1,3));  % asymetric unc.
        end
    end
end
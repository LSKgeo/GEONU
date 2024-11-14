function abundance = abundOut(center,Vp,param,correl_SiO2,correl_abund)
% amph.U.fit.param
% ABUNDOUT will return a value of "X" based on the probability of each X 
% value as defined by "P".  Thus X and P need to be the same length. See 
% this link for more information:
%<a href="matlab:system('start https://www.mathworks.com/matlabcentral/answers/23319-easy-question-with-probability')">https://www.mathworks.com/matlabcentral/answers/23319-easy-question-with-probability</a>
%
% abundOut(probability,center,Vp,param,correl)
%
% correl_sio2 is inputed as a 1x(length) matrix, where "length" is equal to
% the number of SiO2 bins. 
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                  Updated March, 2018                   -----


% - Find index of closest bin for Vp - 
for i = 1:length(center)
    [~,index(i,1)] = min(abs(center(i,:)-Vp(i)));
    P(i,1) = correl_SiO2(i,index(i));  % extract probability of SiO2 for "index" Vp bin
end


% -- Calculate abundance from extracted log-normal parameters -- (exp brings into normal space)
m = exp(param(P,1)); 
eplus = exp(param(P,1)+param(P,2)); % convert log-values to normal
eminus = exp(param(P,1)-param(P,2));

abundance = logdist(m,eplus,eminus,0,correl_abund); % this can be changed to just a MATLAB code
%abundance = exp(randist(param(a,1),param(a,2),0,correl)); 


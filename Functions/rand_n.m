function [output] = rand_n(iterations)
% This function will output "0" if iter = 1 so that we will draw the
% central value of a distribution (instead of randomly along the
% distribution shape). If iter > 1, then we can create random numbers along
% a distribution (normal).  These numbers are single point floating numbers
% because they use less memory. 

    if iterations == 1 
        output = 0;% this will cause the creation of a distribution at the central value
    elseif iterations >1 
        output = randn(iterations,1); 
    end
end
function structure = RM18_allocate(sums,method,cc,layer)

structure.total.hf = stat(sums(:,5),method)/sum(cc); 
structure.total.hp.total = stat(sums(:,6),method); 
structure.total.hp.U = stat(sums(:,7),method); 
structure.total.hp.Th = stat(sums(:,8),method); 
structure.total.hp.K = stat(sums(:,9),method); 
structure.total.mass = stat(sums(:,1),method); 


if strcmp('MC',layer) == 1 | strcmp('LC',layer) == 1 | strcmp('man',layer) == 1
structure.aU = stat(sums(:,2)./(sums(:,1)),method)*10^6;  % aU (ppm)
structure.aTh = stat(sums(:,3)./(sums(:,1)),method)*10^6; % aTh (ppm)
structure.aK = stat(sums(:,4)./(sums(:,1)),method)/0.00011959;  % aK (wt%) 
end





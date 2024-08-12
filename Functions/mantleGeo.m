function [abund_stat_dm,abund_sums_dm,flux_stat_dm,flux_sums_dm,flux_count_dm,...
    abund_stat_em,abund_sums_em,flux_stat_em,flux_sums_em,flux_count_em]...
    = mantleGeo(LAB,lat,lon,aU,aTh,aK40,detector,SurfRadius,PREM,s2)

% MANTLE discretizes the mantle column beneath the 1x1 degree cell defined
% by lat and lon. The code splits the column into smaller voxels depending
% on the distance from the detector.  The equation to voxelize the mantle
% column is the same as the 
%
%
%
% Purpose:  We will take whateve the LAB depth is (central value, ignore
% thickness uncertainty which would make mantle mass uncertainty ~=0.5%).  The mantle goes from LAB to CMB. This mantle we
% will split into ~ 10 layers with thickness of layer increasing with
% depth.  We define the mass of each layer using PREM. We define PREM at
% 1km intervals, allowing us to simply sum up the PREM layers within our
% mantle layer (this accounts for changes in density within our mantle
% layer). 
%
%
%
%  Inputs: 
%  LAB = lithosphere asthenosphere boundary depth (m)
%  lat = latitude of cell center
%  lon = longitude of cell center
%  aU = abundance of Uranium (kg/kg)
%  aTh = abundance of Thorium (kg/kg)
%  aK40 = abundance of Potassium-40 (kg/kg)
%  det = detector location
%  SurfRadius = surface radius (from center of Earth) (m)
%  PREM = Preliminary Reference Earth Model defined every 1km 
%  s2  = Simple2 structure (energies, heat production, calcFlux);


% -- Define Constants --
    CMB_depth = 2891*1e3; % (m) Core Mantle Boundary (CMB) depth

% -- Define Flux parameters -- 
    det = detector(1,:); clear detector% lon2, latidude, SurfRadius (m)

    % - Preallocate variables -
        % for "flux_" variables, NEED to use "zeros" not "nan", otherwise wont sum together
    flux= zeros(1,length(s2.energy));
    temp_flux = zeros(1,19); 
    s2.energy = s2.energy/1000; % convert to MeV for flux calculation
    
    % Define polynomial factors for grid size equation (8th order)
    % The values are ~arbitrary but I found the flux stops changing with these params
    %p = [-2.66418271373068e-42	1.124024966179065e-35	-1.94938406390849e-29	1.78637937516982e-23	-9.23203293588723e-18	2.64355903236367e-12	-3.47340369763152e-07	0.0247642194415582	29.1375158113139]; 

% -- Define bin centers for recording flux vs distance --
    centers = s2.centers; 
    %dist_count_EM = zeros(length(centers),length(s2.energy)); %50x19 matrix

% -- Calculate Area at top of mantle -- (for heatflux calculation)
    % see: http://mathforum.org/library/drmath/view/63767.html
    area_cell = (SurfRadius-LAB).^2.*abs(sind(lat-0.5) - sind(lat+0.5)).*abs((lon-0.5)*pi/180 - (lon+0.5)*pi/180); %(m^2) 

%%


% -- Calculate mass of mantle layers (1 km thick) --
% Note: The masses of layers are calculated assuming layers 1 km thick
% starting from the LAB (which is not an integer). Because I define mantle
% thickness later as integer thicknesses, I can just sum up the masses. 
% The PREM variable is defined from 1:2891, so PREM(100,3) is the density
% at 100 km depth. 
    layers.depth = LAB+500:1000:(CMB_depth-500); %(m) center of each 1km thick layer
    layers.depth = layers.depth'; 
    for i = 1:length(layers.depth)
        layers.mass(i,1) = voxMass(SurfRadius,layers.depth(i),1000,PREM(round(layers.depth(i)/1000+0.5),3),lat-0.5,lat+0.5,lon-0.5,lon+0.5); %kg
        layers.vol(i,1) = layers.mass(i)/PREM(round(layers.depth(i)/1000+500),3); %(m)
    end
    
% -- Calculate mantle thickness --
    thick_man = layers.depth(end)+500 - LAB; %(m)
    
% -- Define thickness of arbitrary mantle layers --
    if s2.calcMantleLayered == true %Enriched and Depleted Mantle
        % Assign arbitrary layer thicknesses. The last layer we will calculate
        % the remaineder after we assign 750km to the Enriched Mantle (~19%
        % mass of mantle (see Arevalo et al. 2013). As 'thick' is defined we
        % will never have a negative last layer (largest LAB depth is 320km). 
        thick = [50,100,150,250,400,400,400,0]; % (km) thickness of layers (arbitrary) 
        thick(end+1) = 750; %=Enriched Mantle Thickness
        thick(end-1) = thick_man/1000 - sum(thick); %thickness of last DM layer = remainder
        thick = thick*1000;

    else %entire mantle = depleted

        % -- Calculate depth and thickness of each layer (~10 layers initial) --
        thick = [50,100,150,250,500,500,500,500]; % (km) thickness of layers (arbitrary) 
        thick(end+1) = thick_man/1000 - sum(thick); %thickness of last layer = remainder
        thick = thick*1000; %(m)
    end


% -- Calculate depth to center of each arbitrary layer --
    depth1 = thick(1)/2+LAB; %first depth
    for i = 2:length(thick)
       depth1(i) = thick(i)/2 + depth1(i-1)+thick(i-1)/2;  % depth to center of layers
    end


% -- Calculate mass (sum masses from layers.mass) -- (kg) 
    x = 1; mass1 = zeros(1,length(thick)); 
    for i = 1:length(thick)
        mass1(i) = sum(layers.mass(x:x+thick(i)/1000-1)); 
        x = x+thick(i)/1000; % move to next layer
    end



% -- Pre-allocate variables --
temp_flux = zeros(1,length(s2.energy),length(mass1)); 
dist_count = zeros(length(centers),length(s2.energy),length(mass1)); %50x19x9 matrix

%temp_flux_EM = zeros(1,length(s2.energy)); 

%% Calculate Heat Production/Flux
    % -- Calculate Mass of HPE in cell (kg) --  
    if s2.calcMantleLayered == true %Enriched and Depleted Mantle
        % Calculate depleted mantle layers (mass1(1:end-1))
        abund_mass(:,1) = aU.dm.*sum(mass1(1:end-1)); % U kg
        abund_mass(:,2) = aTh.dm.*sum(mass1(1:end-1)); % Th kg
        abund_mass(:,3) = aK40.dm.*sum(mass1(1:end-1)); % K40 kg

        % Calculate enriched mantle layers (mass1(end)
        abund_mass(:,4) = aU.em.*sum(mass1(end)); % U kg
        abund_mass(:,5) = aTh.em.*sum(mass1(end)); % Th kg
        abund_mass(:,6) = aK40.em.*sum(mass1(end)); % K40 kg
        
    else 
        % Treat entire mantle as Depleted abundance
        abund_mass(:,1) = aU.dm.*sum(mass1); % U kg
        abund_mass(:,2) = aTh.dm.*sum(mass1); % Th kg
        abund_mass(:,3) = aK40.dm.*sum(mass1); % K40 kg  
        
        % Calculate enriched mantle layers (mass1(end)
        abund_mass(:,4) = zeros(length(aU.dm),1); % U kg
        abund_mass(:,5) = zeros(length(aU.dm),1); % Th kg
        abund_mass(:,6) = zeros(length(aU.dm),1); % K40 kg
         
    end
    
    % -- Calculate Heat Production (W) --
    % - Depleted Mantle 
    heat.U.dm = abund_mass(:,1).* s2.hp.U;
    heat.Th.dm = abund_mass(:,2).* s2.hp.Th;
    heat.K.dm = abund_mass(:,3).*s2.hp.K40;
    heat.total.dm = heat.U.dm + heat.Th.dm + heat.K.dm; 
    
    % - Enriched Mantle 
    heat.U.em = abund_mass(:,4) .* s2.hp.U;
    heat.Th.em = abund_mass(:,5) .* s2.hp.Th;
    heat.K.em = abund_mass(:,6) .*s2.hp.K40;
    heat.total.em = heat.U.em + heat.Th.em + heat.K.em; 
    
    % -- Calculate Heat Flow  (W/m2) --
    heatflow.dm = heat.total.dm./area_cell;
    heatflow.em = heat.total.em./area_cell;


%% Calculate Geoneutrino Flux (only do 1 stage of gridding)

if s2.calcFlux == true 
    for n = 1:length(mass1)
    % ---- Gridding of cells ----
%{
This section seperates the cell into a finer grid equal everywhere within
the cell. For example, if the cell is 100km x 100km x 5 km thick, and I
wish to make cells of 1km wide, then I will produce (100*100*5) = 50,000
new voxels. 

    DO THIS OUTSIDE OF ABOVE LOOP. Otherwise it is too computationally
    intensive.
%}
% -- Calculate distance from center of cell to detector --
%   Remember: this value is sometimes small but the cell size is ~100x100
%   km at equator...

   
%%% GRID LATITUDE AND LONGITUDE %%%

        % - Define gridding size cutoffs - 
            % These cutoffs defines size number of voxels after splitting
    grid.num = [4000,15000,25000,40000,60000]; %[15000,20000,50000,100000,150000]; %(m) ~size of grid cells  [3000,5000,10000,50000,90000]
    grid.lim = [100000,160000,280000,600000,1000000]; %(m) distance from detector for each grid size
     
    distpath = dis(lon,lat,det,SurfRadius,depth1(n)); %(m) 


if distpath > grid.lim(5)
    % - Do not grid cells and calculate geoneutrino flux -
    % We will do all iterations at once
    
    %[mass2,~] = fluxGrid(s1.lon,s1.lat,1,1,depth,thick,rho,SurfRadius,det); 
    
       p1 = 1 - s2.pee.p1 .* sin(1.27.*abs(s2.pee.delm21)*bsxfun(@rdivide,distpath,s2.energy')).^2; 
       p2 = s2.pee.p2 .* sin(1.27.*abs(s2.pee.delm32n).*bsxfun(@rdivide,distpath,s2.energy')).^2;
       p3 = s2.pee.p3 .* sin(1.27.*abs(s2.pee.delm32i).*bsxfun(@rdivide,distpath,s2.energy')).^2;
    
    geoResponse = mass1(n)./(4*pi*distpath.^2); 
    
    unityFlux  = bsxfun(@times,(p1-p2-p3),geoResponse); %kgU/m2 (flux assuming unity concentration)
    
    temp_flux(:,:,n) = unityFlux; 
    
            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath,centers)),[],2);
            dist_count(index,:,n) = dist_count(index,:,n) + unityFlux;     
    
    %{
    if n == length(mass1) && s2.calcMantleLayered == true % loop in Enriched Mantle layer
            temp_flux_EM = unityFlux;
                
            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath,centers)),[],2);
            dist_count_EM(index,:) = dist_count_EM(index,:) + unityFlux; 
        
    else    %treat as depleted mantle
            temp_flux = unityFlux;
            
            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath,centers)),[],2);
            dist_count(index,:) = dist_count(index,:) + unityFlux;             
        
    end
    %}
    
else % Splits the cell into smaller voxels 
    
            % - Calculate lon2 interval -  
    if distpath <= grid.lim(1)
            num = grid.num(1); % num = desired cell size
        elseif distpath > grid.lim(1) && distpath <= grid.lim(2)
            num = grid.num(2);
        elseif distpath > grid.lim(2) && distpath <= grid.lim(3)
            num = grid.num(3);
        elseif distpath > grid.lim(3) && distpath <= grid.lim(4)
            num = grid.num(4);
        elseif distpath > grid.lim(4) && distpath <= grid.lim(5)
            num = grid.num(5); 
    end  
    
       
% Grid latitude and longitude (first gridding)
%(don't do this in loop as lat/lon doesnt change, only depth)


    % Note: WE ARE NOT TAKING INTO ACCOUNT PHYSICAL UNCERTAINTY in this
    % calculation of flux (i.e. we ignore depth and thickness).  We do this
    % because it takes >10 times longer with a negligible change to
    % uncertainty (i think)

     [lon2,lat2,depth2,lon_int2,lat_int2,d_int2]...
        = miniVox(lon,lat,depth1(n),1,1,thick(n),num,SurfRadius); 

    % -- Loop through each new cell and decide if need to grid again --
    %{
      The following code loops through the new cells and based on the
      distance from the detector ('distpath2') will decide if the cell
      needs to be gridded further (see 'num2'). Further gridding will only
      occur very close to the detector
        %}

    for a = 1:length(lon2) 
    
    % - Calc Distance to det for new voxel (m)-                 
    distpath2 = dis(lon2(a),lat2(a),det,SurfRadius,depth2(a)); %(m)
        

            % Calculate Mass of each secondary voxel
        rho = PREM(round(depth2(a)/1000),3); %(kg/m3) Extract PREM densities
        [mass2,~] = fluxGrid(lon2(a),lat2(a),lon_int2,lat_int2,depth2(a),d_int2,rho,SurfRadius,det);     

            % - Calculate Geoneutrino Flux - 
            p1 = 1 - s2.pee.p1 .* sin(1.27.*abs(s2.pee.delm21)*bsxfun(@rdivide,distpath2,s2.energy')).^2; 
            p2 = s2.pee.p2 .* sin(1.27.*abs(s2.pee.delm32n).*bsxfun(@rdivide,distpath2,s2.energy')).^2;
            p3 = s2.pee.p3 .* sin(1.27.*abs(s2.pee.delm32i).*bsxfun(@rdivide,distpath2,s2.energy')).^2;

            geoResponse = mass2./(4*pi*distpath2.^2); 
            unityFlux  = bsxfun(@times,(p1-p2-p3),geoResponse); %kgU/m2 (flux assuming unity concentration)

            
        temp_flux(:,:,n) = temp_flux(:,:,n) + unityFlux; 
        
         % - Record Flux vs Distance (for plotting) -
        [~,index] = min(abs(bsxfun(@minus,distpath2,centers)),[],2);
        dist_count(index,:,n) = dist_count(index,:,n) + unityFlux; 
        
            
            %{
     if n == length(mass1) && s2.calcMantleLayered == true % loop in Enriched Mantle layer
        % - Add up the fluxes from each mini-voxel - 
        temp_flux_EM = temp_flux_EM + unityFlux; 

            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath2,centers)),[],2);
            dist_count_EM(index,:) = dist_count_EM(index,:) + unityFlux; 
    else
        temp_flux = temp_flux  + unityFlux; %dont "sum" as it should be 1 row with xEnergies

            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath2,centers)),[],2);
            dist_count(index,:) = dist_count(index,:) + unityFlux;             
    end
            %}
    end % end of first gridding loop "for" loop

end

% -- Add up flux from each of the mantle layers --
%flux = flux + temp_flux; 
%flux_EM = temp_flux_EM; %only 1 layer so don't need to add like for 'flux'
    end
  
% -- Multiply flux by constants -- This will take the single flux and make
flux_U238.dm = bsxfun(@times,sum(temp_flux(:,:,1:end-1),3),aU.dm); %sum along third dimension
flux_Th232.dm = bsxfun(@times,sum(temp_flux(:,:,1:end-1),3),aTh.dm); 

flux_U238.em = bsxfun(@times,temp_flux(:,:,end),aU.em); % if calcMantleLayered == false, then aU.em = 0
flux_Th232.em = bsxfun(@times,temp_flux(:,:,end),aTh.em); 

% -- Record Flux vs Distance (for plotting) --
dist_count_U238.dm = bsxfun(@times,sum(dist_count(:,:,1:end-1),3),median(aU.dm)); % only do median (otherwise 50x19xiterations)
dist_count_Th232.dm = bsxfun(@times,sum(dist_count(:,:,1:end-1),3),median(aTh.dm));

dist_count_U238.em = bsxfun(@times,dist_count(:,:,end),median(aU.em)); % only do median (otherwise 50x19xiterations)
dist_count_Th232.em = bsxfun(@times,dist_count(:,:,end),median(aTh.em));


%{
% it as long as 'aU'
flux_U238.dm  = bsxfun(@times,flux,aU.dm); 
flux_Th232.dm = bsxfun(@times,flux,aTh.dm); 

flux_U238.em  = bsxfun(@times,flux_EM,aU.em); 
flux_Th232.em = bsxfun(@times,flux_EM,aTh.em); 

% -- Record Flux vs Distance (for plotting) --
dist_count_U238.dm = bsxfun(@times,dist_count,median(aU.dm)); % only do median (otherwise 50x19xiterations)
dist_count_Th232.dm = bsxfun(@times,dist_count,median(aTh.dm)); 

% -- Record Flux vs Distance (for plotting) --
dist_count_U238.em = bsxfun(@times,dist_count_EM,median(aU.dm)); % only do median (otherwise 50x19xiterations)
dist_count_Th232.em = bsxfun(@times,dist_count_EM,median(aTh.dm)); 
%}


end % end of "if geo = 1" 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Organize Data and do Statistics (median + -) for cell 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{ 
    This information is median and +- information for the entire 
    distribution, requiring it to be outside the loop. Calculate statistics
    using "stat" function (created function). 

    I realize putting everything into a single array like "abund_stat" is
    rediculous, but "parfor" doesn't support structures.  This is the best
    way to have a cleaner code/output without having 10 different
    variables. 
    %}
    
    method = 4; %1 = calculate median +-68% c.l.
    % - Mass (kg)
    if s2.calcMantleLayered == true
        abund_stat_dm.mass.mass = stat(sum(mass1(1:end-1)),method); %need 1x3 matrix
        abund_stat_em.mass.mass = stat(mass1(end),method);
    else
        abund_stat_dm.mass.mass = stat(sum(mass1),method);
        abund_stat_em.mass.mass = [0,0,0];  %output is expecting a 1x3 array
    end
    
    % For the below, if calcMantleLayered == false, then the stats for em will be 0
    % - Abundance mass (kg) - 
    abund_stat_dm.mass.U   = stat(abund_mass(:,1),method); 
    abund_stat_dm.mass.Th  = stat(abund_mass(:,2),method);
    abund_stat_dm.mass.K = stat(abund_mass(:,2),method)/.00011959; %kg K
    
    abund_stat_em.mass.U   = stat(abund_mass(:,4),method); 
    abund_stat_em.mass.Th  = stat(abund_mass(:,5),method);
    abund_stat_em.mass.K = stat(abund_mass(:,6),method)/.00011959; %kg K
        
    %- Heat (W)
    abund_stat_dm.hp.U = stat(heat.U.dm,method); 
    abund_stat_dm.hp.Th = stat(heat.Th.dm,method); 
    abund_stat_dm.hp.K = stat(heat.K.dm,method); 
    abund_stat_dm.hp.total = stat(heat.total.dm,method); 
    
    abund_stat_em.hp.U = stat(heat.U.em,method); 
    abund_stat_em.hp.Th = stat(heat.Th.em,method); 
    abund_stat_em.hp.K = stat(heat.K.em,method); 
    abund_stat_em.hp.total = stat(heat.total.em,method); 
    
    %- Heat flow (W/m2)
    abund_stat_dm.hf = stat(heatflow.dm,method); 
    abund_stat_em.hf = stat(heatflow.em,method); 

    

 % --- Combine variables for easy output ---
 %{
     "abund_sums" 
    Col 
     1      mass (kg)
     2      mass of U (kg)
     3      mass of Th (kg) 
     4      mass of K40 (kg)
     5      heat production (W)
     6      heat flow (W/m^2)
 %}
 
     % "abund_sums" are data that undergo iterative summation with every
     % change of cell number (i.e. "n")
     
    if s2.calcMantleLayered == true
        abund_sums_dm(:,1) = repmat(sum(mass1(1:end-1)),length(abund_mass(:,1)),1); %(kg)
        abund_sums_em(:,1) = repmat(mass1(end),length(abund_mass(:,1)),1); %(kg)
    else
        abund_sums_dm(:,1) = repmat(sum(mass1),length(abund_mass(1:3)),1); %(kg)
        abund_sums_em(:,1) = zeros(length(abund_mass(1:3)),1); %(kg)
    end     
    
        % - Depleted Mantle
        abund_sums_dm(:,2) = abund_mass(:,1); % U (kg)
        abund_sums_dm(:,3) = abund_mass(:,2); % Th
        abund_sums_dm(:,4) = abund_mass(:,3); % K40
        abund_sums_dm(:,5) = heatflow.dm; %W/m2
        abund_sums_dm(:,6) = heat.total.dm; % W
        abund_sums_dm(:,7) = heat.U.dm; % W
        abund_sums_dm(:,8) = heat.Th.dm; % W
        abund_sums_dm(:,9) = heat.K.dm; % W
        
        % - Enriched Mantle
        abund_sums_em(:,2) = abund_mass(:,4); % U (kg)
        abund_sums_em(:,3) = abund_mass(:,5); % Th
        abund_sums_em(:,4) = abund_mass(:,6); % K40
        abund_sums_em(:,5) = heatflow.em; %W/m2
        abund_sums_em(:,6) = heat.total.em; % W
        abund_sums_em(:,7) = heat.U.em; % W
        abund_sums_em(:,8) = heat.Th.em; % W
        abund_sums_em(:,9) = heat.K.em; % W        
        
% Preallocate variables LS ALSO CHANGED THIS %%%%%%%%%%%%%%%%%
    %flux_sums_dm = zeros(length(abund_sums_dm),2*length(s2.energy));
    %flux_sums_em = zeros(length(abund_sums_dm),2*length(s2.energy));
    
    flux_stat_dm.U238 = 0; 
    flux_stat_em.U238 = 0; 
    
    flux_stat_dm.Th232 = 0; 
    flux_stat_em.Th232 = 0; 

 
if s2.calcFlux == true
    
   %{
   "flux_sums" are data that undergo iterative summation with every
   change of cell number (i.e. "n") in the main code. 
   
    Col 
   x = length(s2.energy)
     1:x      flux of U238 (neutrino/s)
     x+1:2*x  flux of Th232 (neutrino/s)
   %}
    
   x = length(s2.energy);
   % Geoneutrino Flux
   flux_sums_dm(:,1:x)       = flux_U238.dm; % U238 flux
   flux_sums_dm(:,x+1:2*x)   = flux_Th232.dm;% Th232 flux
   flux_sums_em(:,1:x)       = flux_U238.em; % U238 flux
   flux_sums_em(:,x+1:2*x)   = flux_Th232.em;% Th232 flux   
   
   % Count of flux vs Distance
   flux_count_dm(:,1:x)      = dist_count_U238.dm; 
   flux_count_dm(:,x+1:2*x)  = dist_count_Th232.dm; 
   flux_count_em(:,1:x)      = dist_count_U238.em; 
   flux_count_em(:,x+1:2*x)  = dist_count_Th232.em;      
   
   

   
   %{
   flux_stat is the median +- 68% (confidence limit) for flux of eachisotope for each energy.
   The structure of the variable is first row = median, second row = + uncertainty, 
   and third row is - uncertainty.  Columns are same as "flux_sums"
   
   %}
   % Record Stats about Flux
   for i = 1:length(s2.energy)
       flux_stat_dm.U238  = stat(sum(flux_U238.dm,2),method); 
       flux_stat_em.U238  = stat(sum(flux_U238.em,2),method); 

       flux_stat_dm.Th232 = stat(sum(flux_Th232.dm,2),method);
       flux_stat_em.Th232 = stat(sum(flux_Th232.em,2),method);
   end
end
         
  
function [abund_stat,abund_sums,geoRespStat,flux_stat,flux_sums,temp_pressure,flux_count] = Abund_And_Flux2(k,iter,layer,s1a,s2a,corr,SurfRadius,detector,P,distpath)
%
% ABUND_AND_FLUX calculates the mass, heat production, heat flow, and geoneutrino
% flux for whatever layer you choose. For the middle and lower crust it 
% also calculates the U, Th, and K abundance as well as the other parameters
% listed above using the methodology specified in the "s2" variable. 
%
%   Huang13(k,iter,layer,structure,simplified)
%
%   INPUTS:
% n           = cell number (e.g. cell 3456 out of 64,800)
% iter        = total number of iterations (i.e. distribution length)
% layer       = string of layer name (e.g. 'MC'). Needs to be a string!
% structure   = layer structure (e.g. "UC"), includes 'depth,thick,Vp,pressure,rho,lat,lon'
% s1a          = "simple1" structure containing layer relevant information for specific cell
% s2          = structure containing SurfRadius, hp, oscillation, s2.isotopic, and correlate information.
% cor         = structure with correlation information for layer
% detector    = Detector lat2 and lon2 (in that order)
% P           = Pressure from previous cells combined (length = iter)
%
%   OUTPUTS:  [statistics = median + -] [matrix = distribution of possible values]
% abund_stat = iter x 6 matrix containing info that is iteratively summed
% abund_stat = 64800 long matrix with statistics (median + -) for params
%
%     "abund_sums" 
%     Col 
%      1      mass (kg)
%      2      mass of U (kg)
%      3      mass of Th (kg) 
%      4      mass of K40 (kg)
%      5      heat flow (W/m^2)
%      6      heat production (U+Th+K) (W)
%      7      heat production of U (W)
%      8      heat production of Th (W)
%      9      heat production of K (W)
%
%     "abund_stat"
%     Col
%     1-3     mass (kg; median, + uncertainty, - uncertainty)
%     4-6     mass of U (kg)
%     7-9     mass of Th (kg) 
%     10-12   mass of K40 (kg)
%     13-15   abundance of U (kg/kg)
%     16-18   abundance of Th (kg/kg)
%     19-21   abundance of K40 (kg/kg)
%     22-24   heat production (W)
%     25-27   heat flow (W/m^2)
%     28-30   "f" fraction of felsic (only exists for MC and LC)
%     31-33   temperature in center of layer (only exists for MC and LC)
%     34-36   # of times repeated (only exists for MC and LC)
%
%   REFERENCE for important texts: 
% Huang, Y., Chubakov, V., Mantovani, F., Rudnick, R.L., McDonough, W.F., 2013.
%      A reference Earth model for the heat-producing elements and associated 
%      geoneutrino flux. Geochem. Geophys. Geosystems 14, 2003–2029.
%      <a href="matlab:system('start https://doi.org/10.1002/ggge.20129')">https://doi.org/10.1002/ggge.20129</a>
%
% Dye, S.T., 2012. Geoneutrinos and the radioactive power of the Earth.
%      Reviews of Geophysics 50. 
%      <a href="matlab:system('start https://doi.org/10.1029/2012RG000400')">https://doi.org/10.1029/2012RG000400</a>
%
%
%
%   -----            Written by Scott A. Wipperfurth             ----- 
%   -----      University of Maryland-College Park, Geology      ----- 
%   -----                  Updated January, 2019                 -----



%% Information regarding parfor
%{
temporary variables are those that are used within the parfor but not
used outside of it. For example, I get a temperature value from a
predefined distribution and use it within the parfor iteration, but it is
not necessary to have it outside the parfor environment nor transfer the
specific value between iterations. Thus every iteration re-writes the
value. See https://www.mathworks.com/help/distcomp/temporary-variable.html

%}

%%
% k = 16513;
% layer = 'UC';
% s1a = UC; 
% s2a = simple2;
% corr = cor.UC;
% SurfRadius = GeoPhys.r(k);
% detector = det;
% P = temp_P;
% nf_temp = nearField.logic; % Needed, otherwise "nearField" becomes broadcast variable
% distpath = dis(s1a.lon,s1a.lat,det(1,:),GeoPhys.r,s1a.depth); %(m) s1.depth and det(3) must be depth from surface (not radius)
% distpath = distpath(k);


%% Define constant distributions
 
% -- Check if thickness is 0.  If yes, then dis-continue function --
if s1a.thick(k) == 0  | strcmp('LM_oc',layer) == 1
    % If no layer thickness, then dont do any calculation, just set outputs
    % to 0. Should only apply to sediment layers sometimes. We need to make
    % all of the below variables otherwise the master script will break
    % (due to how structures work)
    
    % Make empty structure.  This is needed so the output structure has
    % same fields as normal (otherwise we get an error)
    abund_stat = zeros(15,3);
    temp_pressure = zeros(iter,1); 
    
    abund_sums = zeros(iter,9);
    geoRespStat = zeros(iter,1);
    
    flux_sums = zeros(iter,2*length(s2a.energy));
    flux_count = zeros(length(s2a.centers),2*length(s2a.energy));
    flux_stat.U238 = 0; 
    %flux_stat.U235 = 0; 
    flux_stat.Th232 = 0; 
    %flux_stat.K40 = 0 ;
    return %this skips the remainder of the function
end



% Convert aK to aK40
s1a.aK40(k,:) = s1a.aK(k,:)*0.00011959; % 0.00011959 converts K to K40 (by mass). Molar is 0.000117



% -- Abundance equation (Huang et al. 2013, eqn 3) --
    abundance = @(aFel,aMaf,f) (aFel .* f)+(aMaf .*(1-f)); %function

    
% -- Set lat/lon bounds of 3D shape (1 degree total boundaries)************
    lonleft = s1a.lon(k) - 0.5;    lonright = s1a.lon(k) + 0.5;
    latbot  = s1a.lat(k) - 0.5;    lattop   = s1a.lat(k) + 0.5;
    r = pi/180;

    
% -- Calculate area_cell of surface above cell -- (for heat flux calculation) --
    % see: http://mathforum.org/library/drmath/view/63767.html
area_cell = SurfRadius.^2.*(sin(latbot*r) - sin(lattop*r)).*(lonleft*r - lonright*r); %(m^2) 
          

% - Preallocate variables -
abund_mass = nan(iter,3); 



% - Extract Bivariate information - 
amp = s2a.bivar.amp; %ok to redefine as it is small size
gran = s2a.bivar.gran; 
seis = s2a.bivar.seis; 


% ----- Geoneutrino Flux Informatoin  ------
% -- Redefine detector information --
if s2a.calcFlux == true
    det = detector(1,:);% clear detector% lon2, latidude, SurfRadius (m)

    % - Preallocate variables -
    geoResponse = zeros(iter,1); 
        % for "flux_" variables, NEED to use "zeros" not "nan", otherwise wont sum together
    flux= zeros(1,length(s2a.energy));
    s2a.energy = s2a.energy/1000; % convert to MeV for flux calculation
    
    % Define polynomial factors for grid size equation (8th order)
    % The values are ~arbitrary but I found the flux stops changing with these params
    p = [-2.66418271373068e-42	1.124024966179065e-35	-1.94938406390849e-29	1.78637937516982e-23	-9.23203293588723e-18	2.64355903236367e-12	-3.47340369763152e-07	0.0247642194415582	29.1375158113139]; 

    
    % -- Define bin centers for recording flux vs distance --
    centers = s2a.centers;
    dist_count = zeros(length(centers),length(s2a.energy)); %50x19 matrix
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% --- Perform matrix calculations ---
%{
    Within this section we will perform calculations based on whether the
    layer is MC, LC, or other.  No matter what layer or type, we will create
    distrbitions of thickness, depth, rho (density), and mass, as these are 
    needed for calculating heat production. These variables are re-written
    every iteration, while variables with  after them are output
    variables. 
    
    For depth, we add 12% uncertainty to the central value. This takes into
    account the uncertainty on the layers above because they all sample the
    distribution at the same location (thanks to the variable "c1".
    %}


% ************* edit the uncertainty**************** LS
    unc = 0.12; % 12% uncertainty on thickness (effects "thick","depth", and "pressure")

% -- Thickness w/ uncer (12 percent) (temperary variable) -- (m)
    temp = strcmp('Crust1',s2a.model); 
    if strcmp('LM',layer) && temp == 1 %LM in CRUST1 or CRUST2
           LAB = randist(175000,75000,0,corr.thick);  
           Moho = randist(s2a.moho(k),s2a.moho(k)*unc,0,corr.thick); 
           thick = LAB - Moho; thick(thick<0) = 0;
    else
        % For other layers or LM in LITHO1.0
            thick = randist(s1a.thick(k),unc*s1a.thick(k),0,corr.thick); thick(thick<0) = 0; %remove negative values if exist
    end
         
% ************************************************************************%    
    
% -- Calculate depth to center of layer (12% uncertainty) -- (m)
    depth = randist(s1a.depth(k), unc*s1a.depth(k), 0, corr.thick); 
    depth(depth<0) = 0; 

% -- Density w/ uncertainty (3 percent) (temporary variables) -- (kg/m3)
        %- value is correlated to Vp
    rho = randist(s1a.rho(k),s1a.rho(k)*0.05,0,corr.vp);   rho(rho<0) = 0; % remove negative if exist

% -- Calculate Mass in each layer -- (kg) 
        % see "help voxMass". Function finds volume of
        % 3d spherical trapezoid then multiplies by density. 
        % SurfRadius is surface radius, depth is to center of layer, thick is thickness of layer
    mass = voxMass(SurfRadius,depth,thick,rho,latbot,lattop,lonleft,lonright,iter);   
    
    
% -- Calculate Pressure from this cell -- (MPa)
    % This value is used for pressure in MC and LC (so we account for
    % uncertainty of parameters). The pressure from this layer gets
    % outputted. 
    temp_pressure = rho .* thick *9.80665 * 10^-6; % MPa
    
     
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Calculate abundance in MC and LC ------- (see Huang et al. 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if strcmp('LC',layer) == 1 || strcmp('MC',layer) == 1 %check if MC or LC

% This is matrix math, not a for-loop.  Matrix math should be a lot faster.

    % -- Vp w/ uncertainty (3 percent) (temporary variable) --  see Huang et al. 2013 and Olugboji et al. 2017 (figure S3A supplement)
    %       - value is correlated to density
    Vp.layer = randist(s1a.Vp(k),0.03*s1a.Vp(k),0,corr.vp); %(km/s) 

    % -- Calculate temperature (for temp correction) -- (equation used by Huang et al. 2013)
        % original equation from Turcotte and Schubert, 2014 equation 4.31
        % general geotherm. k = 3.35 W/m/K, T_0 = 10 C, h_r = 10 km, q_0 =
        % 60 mW/m^2, q_m = 36 mW/m^2 (or 60 mW/m2 *0.6 from relationship observed from Pollack and Chapman, 1977). 
    temperature = 10+71.6.*(1-exp(-depth/1000./10))+10.7.*depth/1000; % (C)

    % -- Calculate pressure (MPa) (for pressure correction) --
    pressure = P + temp_pressure/2; %(kg/m*s2; MPa)
    %pressure = randist(s1.pressure,unc*s1.pressure,0,cor.thick); %(kg/m*s2; MPa)


    
    
    if strcmp('MC',layer) == 1  % Middle Crust---------------------------
        
        if s2a.meth == 1 % Huang et al. 2013 method
    % Amphibolite ("am") felsic ("f") and mafic ("m") (Huang etal. 2013, tbl 5)

            % -- Define mafic and felsic Vp distributions  -- (Huang etal. 2013, tbl 4)
        % - Middle crust (amphibolite)
        Vp.MC.f = randist(6.34,0.16,0,corr.end.vp); %(km/s) (normal dist.)
        Vp.MC.m = randist(6.98,0.20,0,corr.end.vp); %(km/s) (normal dist.)
        
            % -- Temperature + Pressure Correction (subtract 20 C or 600MPa for experimental conditions)--
           % Temperature correct
        Vp.f = Vp.MC.f - (temperature-20)*4*10^-4; % -4 x 10^-4 km/s per C from Huang et al. 2013
        Vp.m = Vp.MC.m - (temperature-20)*4*10^-4; % -4 x 10^-4 km/s per C

            % Pressure correct
        Vp.f = Vp.f + (pressure-600)*2*10^-4; % 2 x 10^-4 km/s per MPa from Huang et al. 2013
        Vp.m = Vp.m + (pressure-600)*2*10^-4; % 2 x 10^-4 km/s per MPa 
        
            % Calculate fraction of felsic ('f')
        f = (Vp.layer - Vp.m)./(Vp.f - Vp.m); % see Huang et al. 2013 equation 2
        x = f; x(x<0)=0; x(x>1)=1; f = x; %no values below 0 or above 1 (need to use x)(COMMENT OUT?)

            % Calculate abundances
        aU   = abundance(s2a.am.f.U,  s2a.am.m.U,  f);  %aU(aU<0,1)=0; %no below 0
        aTh  = abundance(s2a.am.f.Th, s2a.am.m.Th, f);  %aTh(aTh<0,1)=0;
        aK40 = abundance(s2a.am.f.K,  s2a.am.m.K,  f);  %aK40(aK40<0,1)=0;
        
        elseif s2a.meth ==2 % BIVARIATE METHOD middle crust
            
            % Temperature and Pressure correct center Vp 
            center.Vp = bsxfun(@plus, seis.amp.center.Vp,((temperature-20)*-4*10^-4 + (pressure/10^6-600)*2*10^-4)); 
            
            aU = abundOut(center.Vp,Vp.layer,amp.U.fit.param,corr.bivar.sio2,corr.bivar.abund)*10^-6; 
            aTh = abundOut(center.Vp,Vp.layer,amp.Th.fit.param,corr.bivar.sio2,corr.bivar.abund)*10^-6; 
            aK40 = abundOut(center.Vp,Vp.layer,amp.K2O.fit.param,corr.bivar.sio2,corr.bivar.abund)*10^-2*0.000112*0.83; %convert from K2O to K40 
        end   

        
    elseif strcmp('LC',layer) == 1 % Lower Crust---------------------------

        if s2a.meth == 1
        
            %- Lower crust (granulite)
        Vp.LC.f = randist(6.52,0.19,0,corr.end.vp);  %(km/s) (normal dist.)
        Vp.LC.m = randist(7.21,0.20,0,corr.end.vp);  %(km/s) (normal dist.) 
        
        
             % -- Temperature + Pressure Correction (subtract 20 C or 600MPa for experimental conditions)--
           % Pressure correct
        Vp.f = Vp.LC.f - (temperature-20)*4*10^-4; % -4 x 10^-4 km/s per C from Huang et al. 2013
        Vp.m = Vp.LC.m - (temperature-20)*4*10^-4; % -4 x 10^-4 km/s per C

            % Temp correct
        Vp.f = Vp.f + (pressure-600)*2*10^-4; % 2 x 10^-4 km/s per MPa from Huang et al. 2013
        Vp.m = Vp.m + (pressure-600)*2*10^-4; % 2 x 10^-4 km/s per MPa 
        
        
            % Calculate fraction of felsic ('f')
        f = (Vp.layer - Vp.m)./(Vp.f - Vp.m);
        x = f; x(x<0)=0; x(x>1)=1; f = x; %no values below 0 or above 1 (need to use x)(COMMENT OUT?)
          
        % Calculate abundances
        aU   = abundance(s2a.gr.f.U,  s2a.gr.m.U,  f);  %aU(aU<0,1)=0; %no below 0
        aTh  = abundance(s2a.gr.f.Th, s2a.gr.m.Th, f);  %aTh(aTh<0,1)=0;
        aK40 = abundance(s2a.gr.f.K,  s2a.gr.m.K,  f);  %aK40(aK40<0,1)=0;
        
        elseif s2a.meth ==2  % BIVARIATE ANALYSIS
                        % Temperature and Pressure correct center Vp 
            center.Vp = bsxfun(@plus, seis.gran.center.Vp,((temperature-20)*-4*10^-4 + (pressure/10^6-600)*2*10^-4)); 
            
            aU = abundOut(center.Vp,Vp.layer,gran.U.fit.param,corr.bivar.sio2,corr.bivar.abund)*10^-6; 
            aTh = abundOut(center.Vp,Vp.layer,gran.Th.fit.param,corr.bivar.sio2,corr.bivar.abund)*10^-6; 
            aK40 = abundOut(center.Vp,Vp.layer,gran.K2O.fit.param,corr.bivar.sio2,corr.bivar.abund)*10^-2*0.000112*0.83; %convert from K2O to K40 

        end
    end
     
    
else %  Assign abundances if not in "MC" or "LC"
    if s1a.aU(k,2) ==s1a.aU(k,3) % = normal distribution
        aU = randist(s1a.aU(k,1),s1a.aU(k,2),0,corr.abund);
        aTh = randist(s1a.aTh(k,1),s1a.aTh(k,2),0,corr.abund);
        aK40 = randist(s1a.aK40(k,1),s1a.aK40(k,2),0,corr.abund);
        
        aU(aU<0,1)=0; aTh(aTh<0,1)=0; aK40(aK40<0,1)=0; 
        
    else  % = log-normal distribution (only for LM I think?)
        aU = logdist(s1a.aU(k,1),s1a.aU(k,2),s1a.aU(k,3),0,corr.abund);
        aTh = logdist(s1a.aTh(k,1),s1a.aTh(k,2),s1a.aTh(k,3),0,corr.abund);
        aK40 = logdist(s1a.aK40(k,1),s1a.aK40(k,2),s1a.aK40(k,3),0,corr.abund);
        
        aU(aU<0,1)=0; aTh(aTh<0,1)=0; aK40(aK40<0,1)=0; 

    end
   


end %end of loop for calculating abundance in MC or LC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the Th/U ratio and recalculate values if Th/U is outside set bounds

% sigmaratio = 1;
% medratio = 3.77;
% 
% %counter=0;
% 
% %While Th/U value is too high or low, regenerate Th content
% %while
% if aTh/aU > medratio+sigmaratio | aTh/aU < medratio-sigmaratio
%     %counter = counter+1
%     aTh = aU*medratio + (randn(1).*sigmaratio)*(1E-6);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % -- Calculate Mass of HPE in cell (kg) --  
    abund_mass(:,1) = aU  .* mass; % U kg
    abund_mass(:,2) = aTh .* mass; % Th kg
    abund_mass(:,3) = aK40.* mass; % K40 kg


    % -- Calculate Heat Production (W) --
    heat(:,1) = mass .* aU .* s2a.hp.U;
    heat(:,2) = mass .* aTh .* s2a.hp.Th;
    heat(:,3) = mass .* aK40 .*s2a.hp.K40;
    heat(:,4) = sum(heat(:,1:3),2);

    % heat.U = mass .* aU .* s2.hp.U;
    % heat.Th = mass .* aTh .* s2.hp.Th;
    % heat.K = mass .* aK40 .*s2.hp.K40;
    % heat.total = heat.U + heat.Th + heat.K; 
    %{
    if mass <0 % this shouldnt be possible since I stipulate abundances > 0 
        warning('warning heat production <0')
        fprintf('%d %s %d %d %d %d',n,layer,aU,aTh,aK40,mass)
    end
    %}
    % -- Calculate Heat Flow  (W/m2) --
    heatflow = heat(:,4)./area_cell;
%    heatflow = heat.total./area_cell;
      
    

%% Calculate Flux 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Geoneutrino Flux Calculation f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if s2a.calcFlux == true % geo = 1 means calculate flux 
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
%distpath = dis(s1.lon(k),s1.lat(k),det,SurfRadius,s1.depth(k)); %(m) s1.depth and det(3) must be depth from surface (not radius)

%%% GRID LATITUDE AND LONGITUDE %%%

        % - Define gridding size cutoffs - 
            % These cutoffs defines size number of voxels after splitting
    grid.num = [4000,15000,25000,40000,60000]; %[15000,20000,50000,100000,150000]; %(m) ~size of grid cells  [3000,5000,10000,50000,90000]
    grid.lim = [100000,160000,280000,600000,1000000]; %(m) distance from detector for each grid size
    
    
    

if distpath > grid.lim(5)
    % Note: WE ARE NOT TAKING INTO ACCOUNT PHYSICAL UNCERTAINTY in this
    % calculation of flux (i.e. we ignore depth and thickness).  We do this
    % because it takes >10 times longer with a negligible change to
    % uncertainty. Physical uncertainty is indirectly accounted for by the
    % mass of the cell.

    p1 = 1 - s2a.pee.p1 .* sin(1.27.*abs(s2a.pee.delm21)*bsxfun(@rdivide,distpath,s2a.energy')).^2; 
    p2 = s2a.pee.p2 .* sin(1.27.*abs(s2a.pee.delm32n).*bsxfun(@rdivide,distpath,s2a.energy')).^2;
    p3 = s2a.pee.p3 .* sin(1.27.*abs(s2a.pee.delm32i).*bsxfun(@rdivide,distpath,s2a.energy')).^2;

    geoResponse = mass./mass./(4*pi*distpath.^2); % WHAT IS MASS/MASS ABOUT
    flux  = bsxfun(@times,(p1-p2-p3),geoResponse); %kgU/m2 (flux assuming unity concentration)

    % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath,centers)),[],2);
            dist_count(index,:) = dist_count(index,:) + median(flux); 
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
 else % Grid latitude and longitude (first gridding)

            % - Calculate lon2 interval -  
    if distpath <= grid.lim(1)
            num = grid.num(1); % num = desired cell size
        elseif distpath > grid.lim(1) && distpath <= grid.lim(2)
            num = grid.num(2);
        elseif distpath > grid.lim(2) && distpath <= grid.lim(3)
            num = grid.num(3);
        elseif distpath > grid.lim(3) && distpath <= grid.lim(4)
            num = grid.num(4);
        elseif distpath > grid.lim(4) && distpath <= grid.lim(5) %xx
            num = grid.num(5); 
    end  
    
    
    
     [lon2,lat2,depth2,lon_int2,lat_int2,d_int2]...
        = miniVox(s1a.lon(k),s1a.lat(k),s1a.depth(k),1,1,s1a.thick(k),num,SurfRadius); 

    % -- Loop through each new cell and decide if need to grid again --
    %{
      The following code loops through the new cells and based on the
      distance from the detector ('distpath2') will decide if the cell
      needs to be gridded further (see 'num2'). Further gridding will only
      occur very close to the detector
        %}

    for a = 1:length(lon2) 

    % - Calc Distance to det for new small voxel (m)-                 
    distpath2 = dis(lon2(a),lat2(a),det,SurfRadius,depth2(a)); %(m)
    
    % - Calc size of sides of new voxel (m)- 
    num2 = polyval(p,distpath2); %(m) 8th order polynomial with parameters 'p'
    num2(num2<100)=100; % <100 m size doesnt improve signal but takes longer
    
        %Check if desired voxel size is smaller than actual
        if num2*1.333 <= num %i.e. is voxel larger than desired size.  If the new size*1.33 isnt larger than it isnt worth gridding


            % -- Perform a tertiary gridding --
            [lon3,lat3,depth3,lon_int3,lat_int3,d_int3]...
            = miniVox(lon2(a),lat2(a),depth2(a),lon_int2,lat_int2,d_int2,num2,SurfRadius);   

            [mass3,distpath3] = fluxGrid(lon3,lat3,lon_int3,lat_int3,depth3,d_int3,s1a.rho(k),SurfRadius,det); 

            p1 = 1 - s2a.pee.p1 .* sin(1.27.*abs(s2a.pee.delm21)*bsxfun(@rdivide,distpath3,s2a.energy')).^2; 
            p2 = s2a.pee.p2 .* sin(1.27.*abs(s2a.pee.delm32n).*bsxfun(@rdivide,distpath3,s2a.energy')).^2;
            p3 = s2a.pee.p3 .* sin(1.27.*abs(s2a.pee.delm32i).*bsxfun(@rdivide,distpath3,s2a.energy')).^2;

            geoResponse = mean(mass3/mean(mass)./(4*pi*distpath3.^2)); 
            unityFlux  = bsxfun(@times,(p1-p2-p3),geoResponse) ; %kgU/m2 (flux assuming unity concentration)

            %[~,unityFlux,unityFlux_Th232] = geoFlux(mass3,distpath3,s2energy,s2.pee);
            %mass_temp = sum(mass3)+ mass_temp; %used when checking if correct

            % - Sum the fluxes from each mini-voxel - 
            flux = flux  + sum(unityFlux); % 'sum' sums flux for each energy

            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath3,centers)),[],2);
            x = unique(index);
            for i = 1:length(x)
            dist_count(x(i),:) = dist_count(x(i),:) + sum(unityFlux(index==x(i),:)); 
            end
        
        
        else %  Calculate flux using first gridding

            % Calculate Mass of each secondary voxel
            [mass2,~] = fluxGrid(lon2(a),lat2(a),lon_int2,lat_int2,depth2(a),d_int2,s1a.rho(k),SurfRadius,det);     

            %mass_temp = sum(mass2)+ mass_temp;   %used when checking if correct
            % - Calculate Geoneutrino Flux - (calc geoResponse later

            p1 = 1 - s2a.pee.p1 .* sin(1.27.*abs(s2a.pee.delm21)*bsxfun(@rdivide,distpath2,s2a.energy')).^2; 
            p2 = s2a.pee.p2 .* sin(1.27.*abs(s2a.pee.delm32n).*bsxfun(@rdivide,distpath2,s2a.energy')).^2;
            p3 = s2a.pee.p3 .* sin(1.27.*abs(s2a.pee.delm32i).*bsxfun(@rdivide,distpath2,s2a.energy')).^2;

            geoResponse = mean(mass2/mean(mass)./(4*pi*distpath2.^2)); 
            unityFlux  = bsxfun(@times,(p1-p2-p3),geoResponse); %kgU/m2 (flux assuming unity concentration)

            % - Add up the fluxes from each mini-voxel - 
            flux = flux  + unityFlux; %dont "sum" as it should be 1 row with xEnergies

            % - Record Flux vs Distance (for plotting) -
            [~,index] = min(abs(bsxfun(@minus,distpath2,centers)),[],2);
            x = unique(index); %find unique centers
            for i = 1:length(x) %loop through unique centers and sum for each
            dist_count(x(i),:) = dist_count(x(i),:) + sum(unityFlux(index==x(i),:)); 
            end
            


        end % end of second gridding "if" statement              
    end % end of first gridding loop "for" loop
    
    % - Calculate geoResponse for entire original voxel instead of for every mini-vox - 
        %geoResponse = mass/(4*pi*distpath^2); 

    
%     
%        if temp_mass ~= mass %masses should be EXACTLY the same!
%        warning('(no grid) Mass and geoFlux mass not the same!!!%.3e vs %.3e \n',median(mass),median(temp_mass))
%        end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Multiply flux by constants -- This will take the single flux and make
% it as long as 'aU'
flux_U238  = bsxfun(@times,flux,aU.*mass); %nu/m2/kg %.* s2.iso.avgd .* s2.iso.dc.U238 ./s2.iso.amu.U238 .*s2.iso.molar.U238; 
flux_Th232 = bsxfun(@times,flux,aTh.*mass); %.* s2.iso.avgd .* s2.iso.dc.Th232 ./s2.iso.amu.Th232  

% -- Record Flux vs Distance (for plotting) --
dist_count_U238 = bsxfun(@times,dist_count,median(aU)); % only do median (otherwise 50x19xiterations which is too big file size)
dist_count_Th232 = bsxfun(@times,dist_count,median(aTh)); 


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
    ridiculous, but "parfor" doesn't support structures.  This is the best
    way to have a cleaner code/output without having 10 different
    variables. 
    %} 
%NO. CONVERT THIS TO A TABLE NOW
% METHODS FOR EACH SECTION
% METHOD 1:
% METHOD 4: 
% Call stat once and that's it. we need to call stat on
% masstotal, umass, thmass, kmass, hptotal, hpu,hpth,hpk,

    % - f (fraction of felsic; HIGHLY INACCURATE)
    if s2a.meth ~=1
       f = zeros(1,iter);
    end

    
    method1 = 1; %1 = calculate median +-68% c.l.
    method2 = 4; % geometric mean

    % - Mass (kg)
    abund_tot_1 = zeros(iter,9);
    abund_tot_1(:,1) = mass; %total mass for each voxel
    abund_tot_1(:,2:5) = heat;
    abund_tot_1(:,6) = heatflow;
    
    %only for MC or LC
    if strcmp('LC',layer) == 1 || strcmp('MC',layer) == 1 %check if MC or LC
        abund_tot_1(:,7) = pressure;
        abund_tot_1(:,8) = f;
        abund_tot_1(:,9) = temperature;
    end


    abund_tot_2 = zeros(iter,6);
    abund_tot_2(:,1) = abund_mass(:,1); %mass U
    abund_tot_2(:,2) = abund_mass(:,2); %mass Th 
    abund_tot_2(:,3) = abund_mass(:,3); %mass K40
    abund_tot_2(:,4) = aU; %aU concentration
    abund_tot_2(:,5) = aTh; %aTh concentration
    abund_tot_2(:,6) = aK40; %aK40 concentration
    

if strcmp('LC_oc',layer) == 1 || strcmp('MC_oc',layer) == 1
    % Make empty structures.  When we do the MC or LC calc for continental
    % crust it creates these. In order to output to the same structure, we
    % need to make these empty variables, otherwise error. 
    as1= zeros(9,3); 
    as2 = zeros(6,3); 
    
else
    as1 = stat2(abund_tot_1,method1);
    as2 = stat2(abund_tot_2,method2);
    
end     
%     "abund_stat"
%     Col FIX THIS KEY
%     1-3     mass (kg; median, + uncertainty, - uncertainty)
%     4-6     mass of U (kg)
%     7-9     mass of Th (kg) 
%     10-12   mass of K40 (kg)
%     13-15   abundance of U (kg/kg)
%     16-18   abundance of Th (kg/kg)
%     19-21   abundance of K40 (kg/kg)
%     22-24   heat production (W)
%     25-27   heat flow (W/m^2)
%     28-30   "f" fraction of felsic (only exists for MC and LC)
%     31-33   temperature in center of layer (only exists for MC and LC)
%     34-36   # of times repeated (only exists for MC and LC)
abund_stat(1:9,1:3) = as1;
abund_stat(10:15,1:3) = as2;


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
        abund_sums(:,1) = mass; %(kg)
        abund_sums(:,2) = abund_mass(:,1); % U (kg)
        abund_sums(:,3) = abund_mass(:,2); % Th (kg)
        abund_sums(:,4) = abund_mass(:,3); % K40 (kg)
        abund_sums(:,5) = heatflow; %W/m2
        abund_sums(:,6) = heat(:,4);%heat.total; % W
        abund_sums(:,7) = heat(:,1); %heat.U; % W from U
        abund_sums(:,8) = heat(:,2);%heat.Th; % W from Th
        abund_sums(:,9) = heat(:,3);%heat.K; % W from K
        
        
% Preallocate variables
    %geoRespStat = zeros(iter,3);
    
    
    flux_sums = zeros(iter,2*length(s2a.energy));
    flux_stat.U238 = 0; 
    %flux_stat.U235 = 0; 
    flux_stat.Th232 = 0; 
    %flux_stat.K40 = 0 ;
 
if s2a.calcFlux == true
    
   %{
   "flux_sums" are data that undergo iterative summation with every
   change of cell number (i.e. "n")
   
    Col 
   x = length(s2.energy)
     1:x      flux of U238 (neutrino/s)
     x+1:2*x  flux of Th232 (neutrino/s)
   %}
    
   x = length(s2a.energy);
   flux_sums(:,1:x)       = flux_U238; % U238 flux
   flux_sums(:,x+1:2*x)   = flux_Th232;% Th232 flux

   flux_count(:,1:x)      = dist_count_U238; 
   flux_count(:,x+1:2*x)  = dist_count_Th232; 
   
   %{
   flux_stat is the median +- 68% (confidence limit) for flux of eachisotope for each energy.
   The structure of the variable is first row = median, second row = + uncertainty, 
   and third row is - uncertainty.  Columns are same as "flux_sums"
   
   %}
   for i = 1:length(s2a.energy)
       %not calling stat2 here. have memory to save and call stat2 later.
    flux_stat.U238  = sum(flux_U238,2);  
    flux_stat.Th232 = sum(flux_Th232,2);

   end

   % - Record geoResponse value -
   geoRespStat = geoResponse; % in case you had to minivox
   
end
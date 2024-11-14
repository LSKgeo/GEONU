function [simple2, E_int, sigma, weak] = nuOscParam(simple2)

% this is just a simple cut and paste of the original script that loads and
% defines isotope heat production info and geoneutrino oscillation
% parameters. It was put into a function to shorten the main script, but as
% of now requires no inputs.

%    % Isotope information (hp, mass, molar ratio)
% % Define heat Production from kg of isotope (McDonough et al., 2020)
%     simple2.hp.U238 = 94.936*10^-6;  %(Watt/kg)
%     simple2.hp.U235 = 561.65*10^-6; %(Watt/kg)
%     simple2.hp.U    = 98.293*10^-6;  %(Watt/kg)
%     simple2.hp.Th   = 26.180*10^-6;  %(Watt/kg)
%     simple2.hp.K40  = 29.029*10^-6;  %(Watt/kg) K40, not K
% 
% % Define isotope parameters
% % Molar ratios from Chart of Nuclides (National Nuclear Data Center). Isotope masses from Table of Nuclides 
%     % (Nuclear Data Center at KAERI). Avogadro's Number from NIST.
% 
%     % - Molar Ratio - 
%     simple2.iso.molar.U238  = 0.992740; 
%     simple2.iso.molar.U235  = 1 - simple2.iso.molar.U238;
%     simple2.iso.molar.Th232 = 1.00; 
%     simple2.iso.molar.K40   = 0.0001167; % mass ratio = 0.000112
% 
%     % - Isotope Mass - (amu)
%     simple2.iso.amu.U238  = 238.0507826; 
%     simple2.iso.amu.U235  = 235.0439231;
%     simple2.iso.amu.Th232 = 232.0380504; 
%     simple2.iso.amu.K40   = 39.9639987; 
%     simple2.iso.amu.amuKg = 1.660539040*10^-27;
% 
%     % - Avogadro's Number - (atom/mole)
%     simple2.iso.avgd = 6.022140857*10^23; % http://physics.nist.gov/cgi-bin/cuu/Value?na    
% 
%     % - Decay constant - %(/s)
%     simple2.iso.dc.U238  = log(2)/(4.4683e9*3.15569e7);  %(/s)
%     simple2.iso.dc.U235  = log(2)/(0.70348e9*3.15569e7); %(/s) 
%     simple2.iso.dc.Th232 = log(2)/(14.1e9*3.15569e7);   %(/s)
%     simple2.iso.dc.K40   = log(2)/(1.262e9*3.15569e7);   %(/s)

    % Antineutrino parameters 
    % Oscillation, spectrum, cross-section
    % Define geoneutrino oscillation parameters: Values from Tanabashi et al. 2018 (Review of Particle Physics), Patrignani 2016 (Chinese Physics C), Navas et al. 2024 (Particle Data Group), and Capozzi et al. 2017. There is no correlation among parameters (or at least it is not known).
    x = (asind(sqrt(0.0240)) - asind(sqrt(0.0215)))/3; % uncertainty (symmetrical)
    the13               = randist(asind(sqrt(0.0215)),x,1);        %(degree) theta 13 center, +, -
    x = (asind(sqrt(0.354)) - asind(sqrt(0.297)))/3; %positive uncertainty (1 sigma)
    y = (asind(sqrt(0.297)) - asind(sqrt(0.250)))/3; % negative uncertainty
    the12               = logdist(asind(sqrt(0.297)),x,y,1);        %(degree) theta 12 
    simple2.pee.delm21  = logdist(7.37,0.59,0.44,1)*10^-5;  %(eV^2) delta-mass 21 
    simple2.pee.delm32n = logdist(2.56,0.13,0.11,1)*10^-3;  %(eV^2) delta-mass 32 normal hierarchy 
    simple2.pee.delm32i = logdist(-2.54,0.12,0.12,1)*10^-3;  %(eV^2) delta-mass 32 inverted hierarchy 
    weak                = asind(sqrt(0.23122));                     %(degree) weak mixing angle from Canas et al. 2016 Physics Letter B
   
    % Calculate Oscillation probability constants
    % These survival probabilities do not involve distance or energy, so we calculate them now, independent of the Monte Carlo loop.
    simple2.pee.p1 = cosd(the13(1)).^4 .* sind(2.*the12(1)).^2;
    simple2.pee.p2 = sind(2*the13(1)).^2 .* cosd(2.*the12(1)).^2;
    simple2.pee.p3 = sind(2*the13(1)).^2;        
    % Load and bin geoneutrino spectrum - (Enomoto et al. 2007)
    % Then, we define energies at which we will calculate flux. This has to be here so we can pre-define variables.
    E_int =75; %(keV) Energy bin size
    %simple2.energy(:,1) = 1806+E_int/2:E_int:3276-E_int/2; %(keV) 3272 = max but this is easier

    simple2.energy = 0+(E_int/2):E_int:3300; 
    simple2.energy = simple2.energy'; 
    %simple2.energy(:,1) = 0+E_int/2:E_int:3276-E_int/2; %(keV) 3272 = max but this is easier

    %MASTER.energy.bin = simple2.energy; 
    %MASTER.energy.binsize = E_int;

% Enomoto et al. 2007 antineutrino spectrum is then loaded. The original spectrum was binned at 1 keV, 
% so now it must be re-binned into E_int sized energy bins as defined in line 62 of this code.
load('Enomoto2007_AntineutrinoSpectrum.mat' )

        for i = 1:length(simple2.energy)
            top = round(simple2.energy(i)+0.5*E_int);   % top index
            bot = round(simple2.energy(i)-0.5*E_int+1); % bottom index

            simple2.eno.U238(i,1)  = sum(enomoto(bot:top,4));  
            simple2.eno.U235(i,1)  = sum(enomoto(bot:top,5));
            simple2.eno.Th232(i,1) = sum(enomoto(bot:top,6));  
            simple2.eno.K40(i,1)   = sum(enomoto(bot:top,7));  
        end
% Calculate the number of neutrinos per decay chain (for Sramek et al., 2016, flux equation), and then the antineutrino interaction cross-section (see Dye 2012, eqn 15).
% Neutrinos per decay chain
        simple2.eno.num.U238  = sum(enomoto(:,4)); 
        simple2.eno.num.U235  = sum(enomoto(:,5)); 
        simple2.eno.num.Th232 = sum(enomoto(:,6)); 
        simple2.eno.num.K40   = sum(enomoto(:,7)); 

% Calculate antineutrino interaction cross-section
        me = 0.5109989461; %mass electron (MeV)
        delta = 939.5654133 - 938.2720813; %Mass_Neutron - Mass_Proton (MeV)

% Calculate Cross Section (see Dye 2012 equation 11)(note: me^2 replaces 'me' Dye 2012 equation)
    sigma = zeros(1,length(simple2.energy)); 
    e = simple2.energy/1000; 
    idx = e>1.806;
    sigma(idx) = 9.52*((e(idx)-delta).^2).*sqrt(1-(me^2./(e(idx)-delta).^2))*(10^-44); %cm2 interaction cross section as a function of energy

% Define Bin-edges for "Flux vs Distance (m)" plot
        simple2.centers = logspace(2,7.1,100);

end
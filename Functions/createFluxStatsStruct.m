function [tab] = createFluxStatsStruct(flux, final_cut, stat_method1, stat_method2, stat_method3)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% Sediment (TNU) (lognormal or normal, both fit well and give ~same unc.)
    tab.flux.U238(1,:) = stat(flux.sed.sums.cc.U238(final_cut), stat_method1); 
    tab.flux.Th232(1,:) = stat(flux.sed.sums.cc.Th232(final_cut), stat_method1); 


% UC (TNU) (lognormal distribution)
    tab.flux.U238(2,:) = stat(flux.UC.cc.U238(final_cut),stat_method1); 
    tab.flux.Th232(2,:) = stat(flux.UC.cc.Th232(final_cut),stat_method1); 
    tab.flux.total(2,:) = stat(flux.UC.cc.U238(final_cut)+flux.UC.cc.Th232(final_cut),stat_method1); 


% MC (TNU) (lognormal distribution)
    tab.flux.U238(3,:) = stat(flux.MC.cc.U238(final_cut),stat_method2); 
    tab.flux.Th232(3,:) = stat(flux.MC.cc.Th232(final_cut),stat_method2); 
    tab.flux.total(3,:) = stat(flux.MC.cc.U238(final_cut)+flux.MC.cc.Th232(final_cut),stat_method2); 


% LC (TNU) (lognormal distribution)
    tab.flux.U238(4,:) = stat(flux.LC.cc.U238(final_cut),stat_method2); 
    tab.flux.Th232(4,:) = stat(flux.LC.cc.Th232(final_cut),stat_method2); 
    tab.flux.total(4,:) = stat(flux.LC.cc.U238(final_cut)+flux.LC.cc.Th232(final_cut),stat_method2); 


% LM (TNU) (lognormal distribution)
    tab.flux.U238(5,:) = stat(flux.LM.cc.U238(final_cut),stat_method2); 
    tab.flux.Th232(5,:) = stat(flux.LM.cc.Th232(final_cut),stat_method2); 
    tab.flux.total(5,:) = stat(flux.LM.cc.U238(final_cut)+flux.LM.cc.Th232(final_cut),stat_method2); 


% OC Sed (TNU)
    tab.flux.U238(6,:) = stat(flux.sed.sums.oc.U238(final_cut),stat_method1); 
    tab.flux.Th232(6,:) = stat(flux.sed.sums.oc.Th232(final_cut),stat_method1); 
    tab.flux.total(6,:) = stat(flux.sed.sums.oc.U238(final_cut)+flux.sed.sums.oc.Th232(final_cut),stat_method1); 


% OC (TNU) - this was marked as "log normal" but coded as "normal" by the
% original author
    U238_oc = flux.UC.oc.U238 + flux.MC.oc.U238 + flux.LC.oc.U238; 
    Th232_oc= flux.UC.oc.Th232 + flux.MC.oc.Th232 + flux.LC.oc.Th232; 

    tab.flux.U238(7,:) = stat(U238_oc(final_cut),stat_method1); 
    tab.flux.Th232(7,:) = stat(Th232_oc(final_cut),stat_method1); 
    tab.flux.total(7,:) = stat(U238_oc(final_cut) + Th232_oc(final_cut),stat_method1); 

% Bulk CC (Sed + UC + MC + LC) (TNU) (lognormal distribution)
    flux.bulkCC.sums.cc.U238 = flux.UC.cc.U238 + flux.MC.cc.U238 + flux.LC.cc.U238...
            + flux.sed.sums.cc.U238; 
    flux.bulkCC.sums.cc.Th232 = flux.UC.cc.Th232 + flux.MC.cc.Th232 + flux.LC.cc.Th232...
            + flux.sed.sums.cc.Th232; 

    tab.flux.U238(8,:)  = stat(flux.bulkCC.sums.cc.U238(final_cut),stat_method3); 
    tab.flux.Th232(8,:) = stat(flux.bulkCC.sums.cc.Th232(final_cut),stat_method3); 
    tab.flux.total(8,:) = stat(flux.bulkCC.sums.cc.U238(final_cut) + flux.bulkCC.sums.cc.Th232(final_cut),stat_method3); 


% Bulk Crust (Bulk CC + Sed_oc + OC) (TNU) (lognormal distribution)
    U238_crust = flux.bulkCC.sums.cc.U238 + flux.UC.oc.U238 + flux.MC.oc.U238...
        +flux.LC.oc.U238 + flux.sed.sums.oc.U238;
    Th232_crust = flux.bulkCC.sums.cc.Th232 + flux.UC.oc.Th232 + flux.MC.oc.Th232...
        +flux.LC.oc.Th232 + flux.sed.sums.oc.Th232; 
    
    tab.flux.U238(9,:)  = stat(U238_crust(final_cut),stat_method3); 
    tab.flux.Th232(9,:) = stat(Th232_crust(final_cut),stat_method3); 
    tab.flux.total(9,:) = stat(U238_crust(final_cut) + Th232_crust(final_cut),stat_method3); 

    
% Near Field Crust (NFC) (TNU) (lognormal distribution)
    nf_crust.U238 = flux.s1.nf.U238 + flux.s2.nf.U238 + flux.s3.nf.U238...
        + flux.UC.nf.U238 + flux.MC.nf.U238 + flux.LC.nf.U238; 
    nf_crust.Th232 = flux.s1.nf.Th232 + flux.s2.nf.Th232 + flux.s3.nf.Th232...
        + flux.UC.nf.Th232 + flux.MC.nf.Th232 + flux.LC.nf.Th232; 
    
    tab.flux.U238(10,:)  = stat(nf_crust.U238(final_cut),stat_method3); 
    tab.flux.Th232(10,:) = stat(nf_crust.Th232(final_cut),stat_method3); 
    tab.flux.total(10,:) = stat(nf_crust.U238(final_cut) + nf_crust.Th232(final_cut),stat_method3);     
    
    
% Far Field Crust (FFC) (TNU) (lognormal distribution)
    ff_crust.U238 = U238_crust - nf_crust.U238; %total signal - nearField = Far Field
    ff_crust.Th232 = Th232_crust - nf_crust.Th232; 

    tab.flux.U238(11,:)  = stat(ff_crust.U238(final_cut),stat_method3); 
    tab.flux.Th232(11,:) = stat(ff_crust.Th232(final_cut),stat_method3); 
    tab.flux.total(11,:) = stat(ff_crust.U238(final_cut) + ff_crust.Th232(final_cut),stat_method3);
    
    
% Total Lithosphere  (TNU) (lognormal distribution)
    U238_litho = U238_crust + flux.LM.cc.U238; 
    Th232_litho = Th232_crust + flux.LM.cc.Th232; 

    tab.flux.U238(12,:) = stat(U238_litho(final_cut),stat_method3); 
    tab.flux.Th232(12,:) = stat(Th232_litho(final_cut),stat_method3); 
    tab.flux.total(12,:) = stat(U238_litho(final_cut) + Th232_litho(final_cut),stat_method3); 

    
% Depleted Mantle (DM; TNU) (normal distribution (can't do lognormal cause 0 values)    
    tab.flux.U238(13,:) = stat(flux.man.dm.U238(final_cut),stat_method1); 
    tab.flux.Th232(13,:) = stat(flux.man.dm.Th232(final_cut),stat_method1); 
    tab.flux.total(13,:) = stat(flux.man.dm.U238(final_cut)+flux.man.dm.Th232(final_cut),stat_method1); 

    
% Enriched Mantle (EM; TNU)
    tab.flux.U238(14,:) = stat(flux.man.em.U238(final_cut),stat_method1); 
    tab.flux.Th232(14,:) = stat(flux.man.em.Th232(final_cut),stat_method1); 
    tab.flux.total(14,:) = stat(flux.man.em.U238(final_cut)+flux.man.em.Th232(final_cut),stat_method1); 

    
% Total Flux at Detector (TNU) (lognormal distribution)
    flux.total.U238 = U238_litho + flux.man.dm.U238 + flux.man.em.U238; 
    flux.total.Th232 = Th232_litho + flux.man.dm.Th232 + flux.man.em.Th232;
    
    tab.flux.U238(15,:) = stat(flux.total.U238(final_cut),stat_method3); 
    tab.flux.Th232(15,:) = stat(flux.total.Th232(final_cut),stat_method3); 
    tab.flux.total(15,:) = stat(flux.total.U238(final_cut) + flux.total.Th232(final_cut),stat_method3); 
end
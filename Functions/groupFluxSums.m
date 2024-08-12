function [flux] = groupFluxSums(s1_cc_flux_sums,s2_cc_flux_sums,s3_cc_flux_sums,...
    UC_cc_flux_sums, MC_cc_flux_sums,LC_cc_flux_sums,LM_cc_flux_sums,...
    s1_oc_flux_sums,s2_oc_flux_sums,s3_oc_flux_sums,...
    UC_oc_flux_sums,MC_oc_flux_sums,LC_oc_flux_sums,LM_oc_flux_sums,...
    s1_nf_flux_sums,s2_nf_flux_sums,s3_nf_flux_sums,...
    UC_nf_flux_sums,MC_nf_flux_sums,LC_nf_flux_sums,LM_nf_flux_sums,...
    man_flux_sums_em, man_flux_sums_dm, TNU, U238, Th232,...
    s1_distCount_sums, s2_distCount_sums, s3_distCount_sums,...
    UC_distCount_sums, MC_distCount_sums, LC_distCount_sums, LM_distCount_sums,...
    man_distCount_sums_dm, man_distCount_sums_em, s1_nf_flux_sums_cc, s2_nf_flux_sums_cc,s3_nf_flux_sums_cc,...
    UC_nf_flux_sums_cc, MC_nf_flux_sums_cc, LC_nf_flux_sums_cc, LM_nf_flux_sums_cc)


%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Sediment (TNU)
    sed_crust = s1_cc_flux_sums + s2_cc_flux_sums + s3_cc_flux_sums; 
    flux.sed.sums.cc.U238 = sum(bsxfun(@times,sed_crust(:,U238),TNU.U238),2); 
    flux.sed.sums.cc.Th232 = sum(bsxfun(@times,sed_crust(:,Th232),TNU.Th232),2); 

    sed_ocean = s1_oc_flux_sums + s2_oc_flux_sums + s3_oc_flux_sums; 
    flux.sed.sums.oc.U238 = sum(bsxfun(@times,sed_ocean(:,U238),TNU.U238),2); 
    flux.sed.sums.oc.Th232 = sum(bsxfun(@times,sed_ocean(:,Th232),TNU.Th232),2); 

% UC (TNU)
    flux.UC.cc.U238 = sum(bsxfun(@times,UC_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.UC.cc.Th232 = sum(bsxfun(@times,UC_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.UC.oc.U238 = sum(bsxfun(@times,UC_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.UC.oc.Th232 = sum(bsxfun(@times,UC_oc_flux_sums(:,Th232),TNU.Th232),2); 

% MC (TNU)
    flux.MC.cc.U238 = sum(bsxfun(@times,MC_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.MC.cc.Th232 = sum(bsxfun(@times,MC_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.MC.oc.U238 = sum(bsxfun(@times,MC_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.MC.oc.Th232 = sum(bsxfun(@times,MC_oc_flux_sums(:,Th232),TNU.Th232),2); 

% LC (TNU)
    flux.LC.cc.U238 = sum(bsxfun(@times,LC_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.LC.cc.Th232 = sum(bsxfun(@times,LC_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.LC.oc.U238 = sum(bsxfun(@times,LC_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.LC.oc.Th232 = sum(bsxfun(@times,LC_oc_flux_sums(:,Th232),TNU.Th232),2); 
    
% LM (TNU)
    flux.LM.cc.U238 = sum(bsxfun(@times,LM_cc_flux_sums(:,U238),TNU.U238),2); 
    flux.LM.cc.Th232 = sum(bsxfun(@times,LM_cc_flux_sums(:,Th232),TNU.Th232),2); 

    flux.LM.oc.U238 = sum(bsxfun(@times,LM_oc_flux_sums(:,U238),TNU.U238),2); 
    flux.LM.oc.Th232 = sum(bsxfun(@times,LM_oc_flux_sums(:,Th232),TNU.Th232),2); 
   
% Near Field Crust (oc and cc combined)
    flux.s1.nf.U238 = sum(bsxfun(@times,s1_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.s2.nf.U238 = sum(bsxfun(@times,s2_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.s3.nf.U238 = sum(bsxfun(@times,s3_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.UC.nf.U238 = sum(bsxfun(@times,UC_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.MC.nf.U238 = sum(bsxfun(@times,MC_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.LC.nf.U238 = sum(bsxfun(@times,LC_nf_flux_sums(:,U238),TNU.U238),2); 
    flux.LM.nf.U238 = sum(bsxfun(@times,LM_nf_flux_sums(:,U238),TNU.U238),2); 

    %flux.s1.nf.Th232 = sum(bsxfun(@times,s1_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.s1.nf.Th232 = sum(s1_nf_flux_sums(:,Th232).*TNU.Th232,2); 

    flux.s2.nf.Th232 = sum(bsxfun(@times,s2_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.s3.nf.Th232 = sum(bsxfun(@times,s3_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.UC.nf.Th232 = sum(bsxfun(@times,UC_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.MC.nf.Th232 = sum(bsxfun(@times,MC_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.LC.nf.Th232 = sum(bsxfun(@times,LC_nf_flux_sums(:,Th232),TNU.Th232),2); 
    flux.LM.nf.Th232 = sum(bsxfun(@times,LM_nf_flux_sums(:,Th232),TNU.Th232),2); 
    
   
    flux.s1.nf_cc.U238 = sum(bsxfun(@times,s1_nf_flux_sums_cc(:,U238),TNU.U238),2); 
    flux.s2.nf_cc.U238 = sum(bsxfun(@times,s2_nf_flux_sums_cc(:,U238),TNU.U238),2); 
    flux.s3.nf_cc.U238 = sum(bsxfun(@times,s3_nf_flux_sums_cc(:,U238),TNU.U238),2); 
    flux.UC.nf_cc.U238 = sum(bsxfun(@times,UC_nf_flux_sums_cc(:,U238),TNU.U238),2); 
    flux.MC.nf_cc.U238 = sum(bsxfun(@times,MC_nf_flux_sums_cc(:,U238),TNU.U238),2); 
    flux.LC.nf_cc.U238 = sum(bsxfun(@times,LC_nf_flux_sums_cc(:,U238),TNU.U238),2); 
    flux.LM.nf_cc.U238 = sum(bsxfun(@times,LM_nf_flux_sums_cc(:,U238),TNU.U238),2); 

    flux.s1.nf_cc.Th232 = sum(bsxfun(@times,s1_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
    flux.s2.nf_cc.Th232 = sum(bsxfun(@times,s2_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
    flux.s3.nf_cc.Th232 = sum(bsxfun(@times,s3_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
    flux.UC.nf_cc.Th232 = sum(bsxfun(@times,UC_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
    flux.MC.nf_cc.Th232 = sum(bsxfun(@times,MC_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
    flux.LC.nf_cc.Th232 = sum(bsxfun(@times,LC_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
    flux.LM.nf_cc.Th232 = sum(bsxfun(@times,LM_nf_flux_sums_cc(:,Th232),TNU.Th232),2); 
   
    
    
    
% Flux vs Distance Counter
    flux.s1.count.U238 = sum(bsxfun(@times,s1_distCount_sums(:,U238),TNU.U238),2); 
    flux.s2.count.U238 = sum(bsxfun(@times,s2_distCount_sums(:,U238),TNU.U238),2); 
    flux.s3.count.U238 = sum(bsxfun(@times,s3_distCount_sums(:,U238),TNU.U238),2); 
    flux.UC.count.U238 = sum(bsxfun(@times,UC_distCount_sums(:,U238),TNU.U238),2); 
    flux.MC.count.U238 = sum(bsxfun(@times,MC_distCount_sums(:,U238),TNU.U238),2); 
    flux.LC.count.U238 = sum(bsxfun(@times,LC_distCount_sums(:,U238),TNU.U238),2); 
    flux.LM.count.U238 = sum(bsxfun(@times,LM_distCount_sums(:,U238),TNU.U238),2); 
    flux.man.dm.count.U238 = sum(bsxfun(@times,man_distCount_sums_dm(:,U238),TNU.U238),2);
    flux.man.em.count.U238 = sum(bsxfun(@times,man_distCount_sums_em(:,U238),TNU.U238),2);

    flux.s1.count.Th232 = sum(bsxfun(@times,s1_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.s2.count.Th232 = sum(bsxfun(@times,s2_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.s3.count.Th232 = sum(bsxfun(@times,s3_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.UC.count.Th232 = sum(bsxfun(@times,UC_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.MC.count.Th232 = sum(bsxfun(@times,MC_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.LC.count.Th232 = sum(bsxfun(@times,LC_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.LM.count.Th232 = sum(bsxfun(@times,LM_distCount_sums(:,Th232),TNU.Th232),2); 
    flux.man.dm.count.Th232 = sum(bsxfun(@times,man_distCount_sums_dm(:,Th232),TNU.Th232),2);
    flux.man.em.count.Th232 = sum(bsxfun(@times,man_distCount_sums_em(:,Th232),TNU.Th232),2);

% Mantle
    flux.man.dm.U238 = sum(bsxfun(@times,man_flux_sums_dm(:,U238),TNU.U238),2); 
    flux.man.dm.Th232 = sum(bsxfun(@times,man_flux_sums_dm(:,Th232),TNU.Th232),2); 
    
    flux.man.em.U238 = sum(bsxfun(@times,man_flux_sums_em(:,U238),TNU.U238),2); 
    flux.man.em.Th232 = sum(bsxfun(@times,man_flux_sums_em(:,Th232),TNU.Th232),2);
end
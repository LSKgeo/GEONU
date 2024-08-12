function [sums] = groupMassSums(s1_cc_sums,s2_cc_sums,s3_cc_sums,UC_cc_sums, MC_cc_sums,LC_cc_sums,LM_cc_sums,...
    s1_oc_sums,s2_oc_sums,s3_oc_sums,UC_oc_sums,MC_oc_sums,LC_oc_sums,LM_oc_sums,...
    s1_nf_sums,s2_nf_sums,s3_nf_sums,UC_nf_sums,MC_nf_sums,LC_nf_sums,LM_nf_sums,...
    man_sums_em, man_sums_dm,s1_nf_sums_cc,s2_nf_sums_cc, s3_nf_sums_cc, UC_nf_sums_cc, ...
    MC_nf_sums_cc, LC_nf_sums_cc, LM_nf_sums_cc)

% Continental Crust % LAURA might be able to shorten into a for-loop
    sums.s1.cc = s1_cc_sums;
    sums.s2.cc = s2_cc_sums; 
    sums.s3.cc = s3_cc_sums;
    sums.UC.cc = UC_cc_sums; 
    sums.MC.cc = MC_cc_sums; 
    sums.LC.cc = LC_cc_sums; 
    sums.LM.cc = LM_cc_sums; 

% Oceanic Crust
    sums.s1.oc = s1_oc_sums;
    sums.s2.oc = s2_oc_sums; 
    sums.s3.oc = s3_oc_sums;
    sums.UC.oc = UC_oc_sums; 
    sums.MC.oc = MC_oc_sums; 
    sums.LC.oc = LC_oc_sums; 
    sums.LM.oc = LM_oc_sums; 

% Near Field
    sums.s1.nf = s1_nf_sums;
    sums.s2.nf = s2_nf_sums; 
    sums.s3.nf = s3_nf_sums;
    sums.UC.nf = UC_nf_sums; 
    sums.MC.nf = MC_nf_sums; 
    sums.LC.nf = LC_nf_sums; 
    sums.LM.nf = LM_nf_sums; 
    
    
    sums.s1.nf_cc = s1_nf_sums_cc;
    sums.s2.nf_cc = s2_nf_sums_cc; 
    sums.s3.nf_cc = s3_nf_sums_cc;
    sums.UC.nf_cc = UC_nf_sums_cc; 
    sums.MC.nf_cc = MC_nf_sums_cc; 
    sums.LC.nf_cc = LC_nf_sums_cc; 
    sums.LM.nf_cc = LM_nf_sums_cc; 
    
% Mantle
    sums.man.em = man_sums_em; 
    sums.man.dm = man_sums_dm;
end
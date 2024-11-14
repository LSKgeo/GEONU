function [FluxTable] = makeFluxTable(valueFun, TNU,s1_cc_flux_sums,s2_cc_flux_sums,s3_cc_flux_sums,...
    UC_cc_flux_sums,MC_cc_flux_sums,LC_cc_flux_sums,LM_cc_flux_sums,s1_oc_flux_sums,s2_oc_flux_sums,...
    s3_oc_flux_sums,UC_oc_flux_sums,MC_oc_flux_sums,LC_oc_flux_sums,man_flux_sums_dm,man_flux_sums_em)

%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
FluxTable = zeros(11,3);

% U238
FluxTable(1,1) = valueFun(sum(s1_cc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(2,1) = valueFun(sum(s2_cc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(3,1) = valueFun(sum(s3_cc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(4,1) = valueFun(sum(UC_cc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(5,1) = valueFun(sum(MC_cc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(6,1) = valueFun(sum(LC_cc_flux_sums(:,1:44).*TNU.U238,2));

crust_flux_sum = s1_cc_flux_sums(:,1:44) + s2_cc_flux_sums(:,1:44) + s3_cc_flux_sums(:,1:44)...
    + UC_cc_flux_sums(:,1:44) + MC_cc_flux_sums(:,1:44) + LC_cc_flux_sums(:,1:44);
FluxTable(7,1) = valueFun(sum(crust_flux_sum.*TNU.U238,2));

FluxTable(8,1) = valueFun(sum(s1_oc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(9,1) = valueFun(sum(UC_oc_flux_sums(:,1:44).*TNU.U238,2));
FluxTable(10,1) = valueFun(sum(LM_cc_flux_sums(:,1:44).*TNU.U238,2));

total_flux_sum = s1_cc_flux_sums(:,1:44) + s2_cc_flux_sums(:,1:44) +  s3_cc_flux_sums(:,1:44)...
    + UC_cc_flux_sums(:,1:44) + MC_cc_flux_sums(:,1:44) + LC_cc_flux_sums(:,1:44) + LM_cc_flux_sums(:,1:44) + s1_oc_flux_sums(:,1:44)...
    + s2_oc_flux_sums(:,1:44) +  s3_oc_flux_sums(:,1:44) + UC_oc_flux_sums(:,1:44) + MC_oc_flux_sums(:,1:44) + LC_oc_flux_sums(:,1:44);
FluxTable(11,1) = valueFun(sum(total_flux_sum.*TNU.U238,2)); % total lithosphere;

FluxTable(12,1) = valueFun(sum(man_flux_sums_dm(:,1:44).*TNU.U238,2));
FluxTable(13,1) = valueFun(sum(man_flux_sums_em(:,1:44).*TNU.U238,2));
total_bse_flux_sum = total_flux_sum + man_flux_sums_dm(:,1:44) +man_flux_sums_em(:,1:44);

FluxTable(14,1) = valueFun(sum(total_bse_flux_sum.*TNU.U238,2)); % total bulk silicate earth

%--------------------------------------------------%

% Th232
FluxTable(1,2) = valueFun(sum(s1_cc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(2,2) = valueFun(sum(s2_cc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(3,2) = valueFun(sum(s3_cc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(4,2) = valueFun(sum(UC_cc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(5,2) = valueFun(sum(MC_cc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(6,2) = valueFun(sum(LC_cc_flux_sums(:,45:88).*TNU.Th232,2));

crust_flux_sum = s1_cc_flux_sums(:,45:88) + s2_cc_flux_sums(:,45:88) + s3_cc_flux_sums(:,45:88)...
    + UC_cc_flux_sums(:,45:88) + MC_cc_flux_sums(:,45:88) + LC_cc_flux_sums(:,45:88);
FluxTable(7,2) = valueFun(sum(crust_flux_sum.*TNU.Th232,2));

FluxTable(8,2) = valueFun(sum(s1_oc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(9,2) = valueFun(sum(UC_oc_flux_sums(:,45:88).*TNU.Th232,2));
FluxTable(10,2) = valueFun(sum(LM_cc_flux_sums(:,45:88).*TNU.Th232,2));

total_flux_sum = s1_cc_flux_sums(:,45:88) + s2_cc_flux_sums(:,45:88) +  s3_cc_flux_sums(:,45:88)...
    + UC_cc_flux_sums(:,45:88) + MC_cc_flux_sums(:,45:88) + LC_cc_flux_sums(:,45:88) + LM_cc_flux_sums(:,45:88) + s1_oc_flux_sums(:,45:88)...
    + s2_oc_flux_sums(:,45:88) +  s3_oc_flux_sums(:,45:88) + UC_oc_flux_sums(:,45:88) + MC_oc_flux_sums(:,45:88) + LC_oc_flux_sums(:,45:88);
FluxTable(11,2) = valueFun(sum(total_flux_sum.*TNU.Th232,2)); % total lithosphere;


FluxTable(12,2) = valueFun(sum(man_flux_sums_dm(:,45:88).*TNU.Th232,2));
FluxTable(13, 2) = valueFun(sum(man_flux_sums_em(:,45:88).*TNU.Th232,2));
total_bse_flux_sum = total_flux_sum + man_flux_sums_dm(:,45:88) + man_flux_sums_dm(:,45:88);

FluxTable(14,2) = valueFun(sum(total_bse_flux_sum.*TNU.Th232,2)); % total bulk silicate earth

%--------------------------------------------------%

dim = 2; % sum along the second dimension

% Total
FluxTable(1,3) = valueFun(sumFluxes(s1_cc_flux_sums,TNU,dim));
FluxTable(2,3) = valueFun(sumFluxes(s2_cc_flux_sums,TNU,dim));
FluxTable(3,3) = valueFun(sumFluxes(s3_cc_flux_sums,TNU,dim));
FluxTable(4,3) = valueFun(sumFluxes(UC_cc_flux_sums,TNU,dim));
FluxTable(5,3) = valueFun(sumFluxes(MC_cc_flux_sums,TNU,dim));
FluxTable(6,3) = valueFun(sumFluxes(LC_cc_flux_sums,TNU,dim));

crust_flux_sum = s1_cc_flux_sums + s2_cc_flux_sums + s3_cc_flux_sums...
    + UC_cc_flux_sums + MC_cc_flux_sums + LC_cc_flux_sums;
FluxTable(7,3) = valueFun(sumFluxes(crust_flux_sum,TNU,dim));

FluxTable(8,3) = valueFun(sumFluxes(s1_oc_flux_sums,TNU,dim));
FluxTable(9,3) = valueFun(sumFluxes(UC_oc_flux_sums,TNU,dim));
FluxTable(10,3) = valueFun(sumFluxes(LM_cc_flux_sums,TNU,dim));

total_flux_sum = s1_cc_flux_sums + s2_cc_flux_sums +  s3_cc_flux_sums...
    + UC_cc_flux_sums + MC_cc_flux_sums+ LC_cc_flux_sums + LM_cc_flux_sums + s1_oc_flux_sums...
    + s2_oc_flux_sums +  s3_oc_flux_sums + UC_oc_flux_sums + MC_oc_flux_sums + LC_oc_flux_sums;

FluxTable(11,3) = valueFun(sumFluxes(total_flux_sum,TNU,dim)); % total lithosphere;

FluxTable(12,3) = valueFun(sumFluxes(man_flux_sums_dm,TNU,2));
FluxTable(13, 3) = valueFun(sumFluxes(man_flux_sums_em,TNU,2));
total_bse_flux_sum = total_flux_sum + man_flux_sums_dm + man_flux_sums_dm;

FluxTable(14,3) = valueFun(sumFluxes(total_bse_flux_sum,TNU,dim)); % total bulk silicate earth





FluxTable = array2table(FluxTable);
FluxTable.Properties.VariableNames = {'U (TNU)', 'Th (TNU)','Total Flux (TNU)'};

FluxTable.Properties.RowNames = {'Sediment 1','Sediment 2','Sediment 3', 'Upper Crust','Middle Crust',...
    'Lower Crust','Total CC','Oceanic Sediment','Oceanic Crust','Lithospheric Mantle','Total Lithosphere',...
    'Depleted Mantle','Enriched Mantle','Total BSE'};


end
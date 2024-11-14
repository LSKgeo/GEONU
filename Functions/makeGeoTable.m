function [GeoTable] = makeGeoTable(valueFun,s1_cc_sums, s2_cc_sums, s3_cc_sums, UC_cc_sums,...
    MC_cc_sums, LC_cc_sums, LM_cc_sums, s1_oc_sums, s2_oc_sums, s3_oc_sums,UC_oc_sums,...
    MC_oc_sums,LC_oc_sums,man_sums_dm,man_sums_em)

%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

GeoTable = zeros(11,9);
GeoTable(1,:) = valueFun(s1_cc_sums,1);
GeoTable(2,:) = valueFun(s2_cc_sums,1);
GeoTable(3,:) = valueFun(s3_cc_sums,1);
GeoTable(4,:) = valueFun(UC_cc_sums,1);
GeoTable(5,:) = valueFun(MC_cc_sums,1);
GeoTable(6,:) = valueFun(LC_cc_sums,1);
GeoTable(7,:) = valueFun(s1_cc_sums + s2_cc_sums + s3_cc_sums + UC_cc_sums +...
    MC_cc_sums + LC_cc_sums,1); % total continental crust
GeoTable(8,:) = valueFun(s1_oc_sums,1);% + valueFun(s2_oc_sums,1)+valueFun(s3_oc_sums,1);
GeoTable(9,:) = valueFun(UC_oc_sums + MC_oc_sums + LC_oc_sums,1);
GeoTable(10,:) = valueFun(LM_cc_sums,1);% + valueFun(LM_oc_sums,1);
GeoTable(11,:) = valueFun(s1_cc_sums + s2_cc_sums +  s3_cc_sums + UC_cc_sums + MC_cc_sums...
    + LC_cc_sums + LM_cc_sums + s1_oc_sums + s2_oc_sums +  s3_oc_sums...
    + UC_oc_sums + MC_oc_sums + LC_oc_sums ,1); % total lithosphere

% Mantle Values
GeoTable(12,:) = valueFun(man_sums_dm);
GeoTable(13,:) = valueFun(man_sums_em);

% Total Values
GeoTable(14,:) =  valueFun(s1_cc_sums + s2_cc_sums +  s3_cc_sums + UC_cc_sums + MC_cc_sums...
    + LC_cc_sums + LM_cc_sums + s1_oc_sums + s2_oc_sums +  s3_oc_sums...
    + UC_oc_sums + MC_oc_sums + LC_oc_sums+ man_sums_dm + man_sums_em);


GeoTable = array2table(GeoTable);
GeoTable.Properties.VariableNames = {'Total Mass (kg)', 'U (kg)','Th (kg)','K40 (kg)','heat flow (W/m^2)',...
    'Total hp (W)', 'U hp (W)','Th hp (W)','K hp (w)'};

GeoTable.Properties.RowNames = {'Sediment 1','Sediment 2','Sediment 3', 'Upper Crust','Middle Crust',...
    'Lower Crust','Total CC','Oceanic Sediment','Oceanic Crust','Lithospheric Mantle','Total Lithosphere',...
    'Depleted Mantle','Enriched Mantle','Total BSE'};

end
function x = sumFluxes(layer, TNU,dim)

    x = sum(layer(:,1:44).*TNU.U238 + layer(:,45:88).*TNU.Th232,dim);
end


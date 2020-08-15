function save_sdppoly(polycell, polyname, vararray, filename)

    for l = 1:length(polycell)
        poly = polycell{l}.(polyname);
        for i = 1:length(poly)
            [coefs, monobase] = coefficients(poly(i),vararray);
            monodeg = [];
            for j = 1:length(coefs)
                monodeg(j,:) = degree(monobase(j),vararray);
            end
            poly_comp{l}{i}.coefs = coefs;
            poly_comp{l}{i}.monodeg = monodeg;
        end
    end
       
    save(filename,'poly_comp')

end
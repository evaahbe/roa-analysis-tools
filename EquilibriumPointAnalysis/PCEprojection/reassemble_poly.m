function fullpoly = reassemble_poly(vararray, poly_comp,l)
    
      
    for i = 1:length(poly_comp{l})
        poly = 0;
        coefs = poly_comp{l}{i}.coefs;
        monodeg = poly_comp{l}{i}.monodeg;

        for j = 1:length(coefs)
            term = 1;
            for k = 1:length(vararray)
                term = term * vararray(k)^monodeg(j,k);
            end
            poly = poly + coefs(j) * term;
        end
        fullpoly(i) = poly;
    end
    fullpoly = fullpoly(:); %return a column vector;
    

end
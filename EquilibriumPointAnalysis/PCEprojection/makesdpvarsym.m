function fxi_sym = makesdpvarsym(fxi,xi,xs)
  

    for j = 1:length(fxi)
        [x1,vx1] = coefficients(fxi(j),[xi]);
        fxi_sym_pre = 0;
        fxi_sym_vec = sym(zeros(1,length(vx1)));
        for g = 1:length(vx1)
            if value(vx1(g)) ==1
                fxi_sym_pre = value(x1(g));
            else
                multiterm = 1;
                for k = 1:length(xi)
                    pow = exponents(vx1(g),xi(k));
                    multiterm = multiterm * xs(k)^pow;
                end          
                fxi_sym_vec(g) = x1(g)*multiterm;
            end
        end
        fxi_sym(j) = fxi_sym_pre+sum(fxi_sym_vec);                      
    end
    fxi_sym = fxi_sym.';


end
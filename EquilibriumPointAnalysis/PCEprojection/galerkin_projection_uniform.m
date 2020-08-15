function coefficients = galerkin_projection_uniform(a,b,order)

    Pc{1} = [1];
    Pc{2} = [1 0];
    Pc{3} = 0.5 * [3 0 -1];
    Pc{4} = 0.5 * [5 0 -3 0];
    Pc{5} = 1/8 * [35 0 -30 0 3];
    Pc{6} = 1/8 * [63 0 -70 0 15 0];
    Pc{7} = 1/16 * [231 0 -315 0 105 0 -5];
    Pc{8} = 1/16 * [429 0 -693 0 315 0 -35 0];
    Pc{9} = 1/128 * [6435 0 -12012 0 6930 0 -1260 0 35];
    Pc{10} = 1/128 * [12155 0 -25740 0 18018 0 -4620 0 315 0];
    Pc{11} = 1/256 * [46189 0 -109395 0 90090 0 -30030 0 3465 0 -63];
    Pc{12} = 1/256 * [88179 0 -230945 0 218790 0 -90090 0 15015 0 -693 0];
    Pc{13} = 1/1024 * [676039 0 -1939938 0 2078505 0 -1021020 0 225225 0 -18018 0 231];
    
    unif_var = [(b-a)/2 (b+a)/2]; %see either my notes or ''﻿Polynomial chaos expansions and stochastic finite-element methods'', Sudret 2014
    
    for i =1:order+1
        coefficients_int = polyint(conv(unif_var,Pc{i}));
        coefficients(i) = diff(polyval(coefficients_int*0.5,[-1,1]))/legendreProduct(i-1,i-1);
    end


end
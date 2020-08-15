function coefficients = galerkin_projection_normal(mean,var,order)

% for the gaussian variable equation see either my notes or ''﻿Polynomial chaos expansions and stochastic finite-element methods'', Sudret 2014
    
    for i =1:order+1
        int_fun = @(z) (mean + var.*z).*Pc_poly(z,i).*1/sqrt(2*pi).*exp(-z.^2/2);
        coefficients(i) = integral(int_fun, -inf,inf)/factorial(i-1);
    end


end

function [Pcpoly] = Pc_poly(z,i)

    Pc{1} = 1;
    Pc{2} = z;
    Pc{3} = z.^2 -1;
    Pc{4} = z.^3 -3.*z;
    Pc{5} = z.^4 -6.*z.^2 +3;
    Pc{6} = z.^5 -10.*z.^3 +15.*z;
    Pc{7} = z.^6 -15.*z.^4 +45.*z.^2 -15;
    Pc{8} = z.^7 -21.*z.^5 +105.*z.^3 -105.*z;
    Pc{9} = z.^8 -28.*z.^6 +210.*z.^4 -420.*z.^2 +105;
    Pc{10} = z.^9 -36.*z.^7 +378.*z.^5 -1260.*z.^3 +945.*z;
    Pc{11} = z.^10 -45.*z.^8 +630.*z.^6 -3150.*z.^4 +4725.*z.^2 -945;
    
    Pcpoly = Pc{i};
end
%function computing products and projections of Hermite basis functions
%with up to a product of nargin-1 basis functions


function hermite_integral = hermiteProduct(varargin)
    
    %the first varargin contains the (LHS) projection poly order, and the rest
    %the polys building the (RHS) product.

    %when referring to the polynomial orders, always add +1 because matlab
    %starts at 1, while orders at 0.
    

    %compute integral
    hermite_fun =@(z) hermite_product(z,varargin);
    
    hermite_integral= integral(hermite_fun, -inf,inf);
    
 
    
end

function [herm_prod] = hermite_product(z,varin)

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
    
    herm_prod = Pc{varin{2}+1}.*Pc{varin{1}+1};
    
    if length(varin) >= 3       
        for i =3:length(varin) %the first one (i=1) is the argument containing the projection index!
            herm_prod = herm_prod.*Pc{varin{i}+1};
        end
    end
    
    herm_prod = herm_prod.*(1/(sqrt(2*pi))).*exp(-(z.^2)/2);
end
%function computing products and projections of Legendre basis functions
%with up to a product of nargin-1 basis functions


function Legendre_integral = legendreProduct(varargin)
    
    %the first varargin contains the (LHS) projection poly order, and the rest
    %the polys building the (RHS) product.

    %when referring to the polynomial orders, always add +1 because matlab
    %starts at 1, while orders at 0.
        

    %define Legendre polynomials
    
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

    
    %build product
    if nargin == 2
        pre_product = Pc{varargin{2}+1};
        
    elseif nargin == 3
        pre_product = conv(Pc{varargin{2}+1},Pc{varargin{3}+1});
    
    elseif nargin > 3
        pre_product = conv(Pc{varargin{2}+1},Pc{varargin{3}+1});
        for i =4:length(varargin) %the first one (i=1) is the argument containing the projection index!
            poly_index = varargin{i}+1; 
            pre_product = conv(pre_product,Pc{poly_index});
        end
    end
    
    %compute integral
    Legendre_open_integral = polyint(conv(pre_product,Pc{varargin{1}+1}));
    
    Legendre_integral = diff(polyval(Legendre_open_integral*0.5,[-1 1]));
    
   
    
end

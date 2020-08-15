%SCALED
%3rd order Taylor (4 terms)

function [ftay] = uncertain_poly_dynamics_TEMPLATE(sys,onHP,xloc,uncer_real)

    % fill in the appropriate uncertain dynamics in poly approx if not poly
    % already
    % insert (mu(i)+uncer_real(i)) for every uncertain parameter!
    
    mu = sys.mu
    x1 = xloc(1);
    x2 = xloc(2);
    x3 = xloc(3);
    xorb1 = onHP.xorb(1);
    xorb2 = onHP.xorb(2);
    xorb3 = onHP.xorb(3);
    k1 = onHP.K(1);
    k2 = onHP.K(2);
    Zi1_1 = onHP.Z(1,1);
    Zi1_2 = onHP.Z(1,2);
    Zi1_3 = onHP.Z(1,3);
    Zi2_1 = onHP.Z(2,1);
    Zi2_2 = onHP.Z(2,2);
    Zi2_3 = onHP.Z(2,3);
    unom = onHP.uorb;
    
    x1dot = 1;
    x2dot = (mu+uncer_real) * 1;
    
    fx = [x1dot;x2dot];
end
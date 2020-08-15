

function [fxpoly] = poly_dynamics_TEMPLATE(x_on_lc,u_on_lc,sys,x,K,Z)
    
    % fill in the approx poly dynamics (if not poly already)
    
    mu = sys.mu;
    x1dot = 1;
    x2dot = mu * 1;
    
    fxpoly = [x1dot;x2dot];
end
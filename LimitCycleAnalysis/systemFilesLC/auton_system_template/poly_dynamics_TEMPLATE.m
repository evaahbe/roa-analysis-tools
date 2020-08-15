

function [fxpoly] = poly_dynamics_TEMPLATE(x,u, sys)
    
    % fill in the approx poly dynamics (if not poly already)
    
    mu = sys.mu;
    x1dot = 1;
    x2dot = mu * 1;
    
    fxpoly = [x1dot;x2dot];
end
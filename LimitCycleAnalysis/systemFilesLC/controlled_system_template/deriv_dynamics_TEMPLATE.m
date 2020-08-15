

function ffx = deriv_dynamics_TEMPLATE(x,u,sys)
    
    mu = sys.mu;
    
    % fill in the dynamic derivative of the true dynamics
    
    x1dotdot = 1;
    x2dotdot = mu * 1;
    
    ffx = [x1dotdot;x2dotdot];
end
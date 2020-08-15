function fx = dynamics_TEMPLATE(x,u,sys)
    
    % fill in the correct true dynamics
    
    mu = sys.mu;
    x1dot = 1;
    x2dot = mu * 1;
    
    fx = [x1dot;x2dot];
end
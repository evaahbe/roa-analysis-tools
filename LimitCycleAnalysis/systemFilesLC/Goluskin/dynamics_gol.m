function fx = dynamics_gol(x,sys)
    
    mu = sys.mu;
    x1dot = -(x(2) - mu*x(1)*(x(1)^2-1/4)*(x(1)^2-1));
    x2dot = x(1);
    
    fx = [x1dot;x2dot];
end
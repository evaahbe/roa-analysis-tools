function fx = dynamics_quadplan(x,sys)
    
    mu = sys.mu;
    x1dot = x(2)*(x(1) + x(2) + mu) + (x(1)^2 + x(2)^2 -1)*(1-x(1));
    x2dot =  -x(1) *(x(1)+ x(2) + 3);
    
    fx = [x1dot;x2dot];
end
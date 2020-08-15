function fx = uncertain_dynamics_vdp(sys,onHP,xloc,uncer_real)
    
    x(1) = xloc(1);
    x(2) = xloc(2);
    mu = sys.mu;
    
    x1dot = x(2);
    x2dot = (mu+uncer_real) * (1 - x(1)^2) * x(2) - x(1);
    
    fx = [x1dot;x2dot];
end
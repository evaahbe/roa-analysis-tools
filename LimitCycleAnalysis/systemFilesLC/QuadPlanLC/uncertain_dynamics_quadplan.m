function fx = uncertain_dynamics_quadplan(sys,onHP,xloc,uncer_real)
    
    x(1) = xloc(1);
    x(2) = xloc(2);
    mu = sys.mu;
    
    x1dot = x(2)*(x(1) + x(2) + (mu+uncer_real)) + (x(1)^2 + x(2)^2 -1)*(1-x(1));
    x2dot =  -x(1) *(x(1)+ x(2) + 3);
    
    
    fx = [x1dot;x2dot];
end
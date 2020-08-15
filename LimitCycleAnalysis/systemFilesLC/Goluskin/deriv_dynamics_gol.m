function ffx = deriv_dynamics_gol(x,sys)
    
    mu = sys.mu;
    x1dot = -x(2) + mu*x(1)*(x(1)^2-1/4)*(x(1)^2-1);
    x2dot = x(1);
    
    x1dotdot = -x2dot + mu*x1dot*(5*x(1)^4 - 3.75*x(1)^2 + 0.25);
    x2dotdot = x1dot;
    
    ffx = [x1dotdot;x2dotdot];
end
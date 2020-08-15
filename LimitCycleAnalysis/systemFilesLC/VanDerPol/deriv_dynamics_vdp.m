function ffx = deriv_dynamics_vdp(x,sys)
    
    mu = sys.mu;
    x1dot = x(2);
    x2dot = mu * (1 - x(1)^2) * x(2) - x(1);
    
    x1dotdot = x2dot;
    x2dotdot = mu * (x2dot - (x2dot*x(1)^2 + x(2)*2*x(1)*x1dot)) - x1dot;
    
    ffx = [x1dotdot;x2dotdot];
end
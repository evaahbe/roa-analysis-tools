function fx = dynamics_vdp_ep(x,mu)
    
    x1dot = x(2);
    x2dot = mu(1) * (1 - x(1)^2) * x(2) - x(1);
    
    fx = -[x1dot;x2dot];
end
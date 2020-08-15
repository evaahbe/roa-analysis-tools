function fx = dynamics_uc(x,mu)

    x1dot = -mu*x(2)+x(1)*(1-x(1)^2-x(2)^2);
    x2dot = x(1)+x(2)*(1-x(1)^2-x(2)^2);
    
    fx = [x1dot;x2dot];
end
function ffx = dynDeriv_uc(x,mu)

    x1dot = -mu*x(2)+x(1)*(1-x(1)^2-x(2)^2);
    x2dot = x(1)+x(2)*(1-x(1)^2-x(2)^2);
    
    x1dotdot = -mu*x2dot +x1dot*(1-x(1)^2-x(2)^2) + x(1)*(-2*x(1)*x1dot-2*x(2)*x2dot);
    x2dotdot = x1dot +x2dot*(1-x(1)^2-x(2)^2) + x(2)*(-2*x(1)*x1dot-2*x(2)*x2dot);
    
    ffx = [x1dotdot;x2dotdot];
end
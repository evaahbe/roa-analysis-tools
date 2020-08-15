function ffx = deriv_dynamics_quadplan(x,sys)
    
    mu = sys.mu;
    x1dot = x(2)*(x(1) + x(2) + mu) + (x(1)^2 + x(2)^2 -1)*(1-x(1));
    x2dot =  -x(1) *(x(1)+ x(2) + 3);
    
    x1dotdot = x2dot*x(1)+x(2)*x1dot + (2*x(2)*x2dot + x2dot*mu + 2*x(1)*x1dot + 2*x(2)*x2dot)*(1-x(1))-x1dot*(x(1)^2 + x(2)^2 -1);
    x2dotdot = -2*x(1)*x1dot -x2dot*x(1)-x(2)*x1dot-3*x1dot;
    
    ffx = [x1dotdot;x2dotdot];
end
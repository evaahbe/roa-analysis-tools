function fx = dynamics_iannelli_INT(t,x,mu)
         
    x1dot = -x(2) - 3/2*x(1)^2 - 1/2*x(1)^3 + mu;
    x2dot = 3*x(1) - x(2) - x(2)^2;

    fx = [x1dot;x2dot];
end
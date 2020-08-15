function fx = dynamics_chesi_INT(t,x,mu,u)
         
        x1dot =  -x(1) + x(2) -mu*(x(1)^2+x(2)^3) + x(1)*u;          
        x2dot = -2*x(2)-x(1)^2 + u;

    fx = [x1dot;x2dot];
end
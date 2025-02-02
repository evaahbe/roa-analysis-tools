
function f = dynamics_ShortPeriod(x,mu,u)


    x1dot = mu(1)*(-1.492*x(1)^3 + 4.239*x(1)^2 + 3.063*10^(-3)*x(1)*x(2) + 6.226*10^(-3)*x(2)^2) -3.236*x(1) + 9.227*10^(-1)*x(2) + (-3.166*10^(-1) +  2.402*10^(-1)*x(1))*u;
    x2dot = mu(2)*(-7.228*x(1)^3 + 18.36*x(1)^2 + 1.103354559437335*x(2)^3) - 45.34*x(1) - 4.372*x(2) + (- 59.99 + 41.5*x(1))*u ; 

    f = [x1dot; x2dot];

end
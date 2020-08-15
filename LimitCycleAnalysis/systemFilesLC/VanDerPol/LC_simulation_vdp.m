clear all


T   = 50;
x0  = [1.706;-0.7749];
mu  = 1;

%options = odeset('RelTol',1e-11,'AbsTol',1e-12);
options = [];

[tvec, xvec] = ode45(@(t,x) simu_dynamics_vdp(t,x,mu),[0 T],x0,options);

figure
plot(xvec(:,1),xvec(:,2))


function fx = simu_dynamics_vdp(t,x,mu)

    x1dot = x(2);
    x2dot = mu * (1 - x(1)^2) * x(2) - x(1);
    
    fx = [x1dot;x2dot];
end




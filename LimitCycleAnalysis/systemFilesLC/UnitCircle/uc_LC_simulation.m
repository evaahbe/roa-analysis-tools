clear all

T =10;

x0 = [-0.915;-0.386];

mu = 1;

%options = odeset('RelTol',1e-11,'AbsTol',1e-12);
options = [];

[tvec, xvec] = ode45(@(t,x) simu_dynamics_vdp(t,x,mu),[0 T],x0,options;

figure
plot(xvec(:,1),xvec(:,2))


function fx = simu_dynamics_vdp(t,x,mu)

    x1dot = -mu*x(2)+x(1)*(1-x(1)^2-x(2)^2);
    x2dot = x(1)+x(2)*(1-x(1)^2-x(2)^2);
    
    fx = [x1dot;x2dot];
end


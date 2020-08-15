% this file is required for simulation of the closed loop system

function f = dynamics_TEMPLATE_INT(t,x,mu1,mu2,K)

    u = K*[x(1),x(2),x(2)^2,...]'; % enter here the chosen control law u = K*h(x)
    x1dot = ;
    x2dot = ; 

    f = [x1dot; x2dot];

end
function [A] = linear_dynamics_vdp(sys, t)
     
    x = ppval(sys.xorb,t);  
    mu = sys.mu;
 
    A = [0, 1; -2 *mu * x(1) * x(2) - 1, mu*(1 - x(1)^2)];
    
end
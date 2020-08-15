function [A] = linear_dynamics_gol(sys, t)
     
    x = ppval(sys.xorb,t);  
    mu = sys.mu;
    
    A = [mu*(5*x(1)^4 - 3.75*x(1)^2 +0.25), -1;...
        1,0];
 
    
end
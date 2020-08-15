function [A] = linear_dynamics_TEMPLATE(sys, t)
     
    x = ppval(sys.xorb,t);  
    mu = sys.mu;
    
    % fill in the correct jacobian
    A = [1,1;1,1]; 
    
end
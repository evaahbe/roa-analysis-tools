function [A] = linear_dynamics_quadplan(sys, t)
     
    x = ppval(sys.xorb,t);  
    mu = sys.mu;
 
%     x1dot = x(2)*(x(1) + x(2) + mu) + (x(1)^2 + x(2)^2 -1)*(1-x(1));
%     x2dot =  -x(1) *(x(1)+ x(2) + 3);
    
    A = [x(2)-3*x(1)^2-x(2)^2+1+2*x(1),2*x(2)+mu+x(1)+2*x(2)-2*x(2)*x(1);-2*x(1)-x(2)-3,-x(1)];
    
end
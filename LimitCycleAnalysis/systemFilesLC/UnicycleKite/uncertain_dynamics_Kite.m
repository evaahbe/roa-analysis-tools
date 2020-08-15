%Scaled


function [fx] = uncertain_dynamics_Kite(x, u, sys, uncer_real)
    
    f1 = sys.paras.vktmax * cos(x(1)) * cos(x(2)) * cos(10*x(3));
    f2 = sys.paras.vktmax * cos(x(2)) * sin(10*x(3));
    f3 = (sys.ugain + uncer_real) * cos(x(1)) * cos(x(2)) * u;
    
    fx = [f1;f2;f3];
    
    
end
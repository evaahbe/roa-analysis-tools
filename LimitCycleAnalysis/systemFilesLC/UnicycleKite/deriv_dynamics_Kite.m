%SCALED

function [dfx] = deriv_dynamics_Kite(x,u,sys,t)
    
    udot = ppval(sys.uorbdot,t);

     
    f1 = sys.paras.vktmax * cos(x(1)) * cos(x(2)) * cos(10*x(3));
    f2 = sys.paras.vktmax * cos(x(2)) * sin(10*x(3));
    f3 = sys.ugain * cos(x(1)) * cos(x(2)) * u;
    

    dfx = [ sys.paras.vktmax  * (-sin(x(1)) * f1 * cos(x(2))* cos(10*x(3)) - sin(x(2)) * f2 * cos(x(1))* cos(10*x(3)) - cos(x(2)) * cos(x(1))* sin(10*x(3))*10*f3);...
           sys.paras.vktmax * (-sin(10*x(3)) * sin(x(2)) * f2 + cos(10*x(3)) *10* f3 * cos(x(2)));...
           sys.ugain * (-sin(x(1)) * f1 * cos(x(2)) * u - sin(x(2)) * f2 * cos(x(1)) * u + cos(x(2)) * cos(x(1)) * udot)];



end
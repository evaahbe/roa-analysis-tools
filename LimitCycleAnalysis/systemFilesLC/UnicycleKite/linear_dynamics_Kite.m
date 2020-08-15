%SCALED

function [A,B] = linear_dynamics_Kite(sys,t)

    x = ppval(sys.xorb,t);
    u = ppval(sys.uorb,t);
    
    A = [-sys.paras.vktmax *sin(x(1)) * cos(x(2)) * cos(10*x(3)), -sys.paras.vktmax * cos(x(1)) * sin(x(2)) * cos(10*x(3)), -sys.paras.vktmax * cos(x(1)) * cos(x(2)) * 10*sin(10*x(3));...
        0, -sys.paras.vktmax * sin(x(2)) * sin(10*x(3)), sys.paras.vktmax * cos(x(2)) * 10*cos(10*x(3));
        -sys.ugain * sin(x(1)) * cos(x(2)) * u, -sys.ugain * cos(x(1)) * sin(x(2)) * u, 0];
    
    B = [0;...
        0;...
        sys.ugain * cos(x(1)) * cos(x(2))];


end
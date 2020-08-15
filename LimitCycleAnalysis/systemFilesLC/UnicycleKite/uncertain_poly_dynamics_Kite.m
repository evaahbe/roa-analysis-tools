%SCALED
%3rd order Taylor (4 terms)

function [ftay] = uncertain_poly_dynamics_Kite(sys,onHP,xloc,uncer_real)
    
    x1 = xloc(1);
    x2 = xloc(2);
    x3 = xloc(3);
    xorb1 = onHP.xorb(1);
    xorb2 = onHP.xorb(2);
    xorb3 = onHP.xorb(3);
    k1 = onHP.K(1);
    k2 = onHP.K(2);
    Zi1_1 = onHP.Z(1,1);
    Zi1_2 = onHP.Z(1,2);
    Zi1_3 = onHP.Z(1,3);
    Zi2_1 = onHP.Z(2,1);
    Zi2_2 = onHP.Z(2,2);
    Zi2_3 = onHP.Z(2,3);
    unom = onHP.uorb;
    
    vktmax = sys.paras.vktmax;
    ugain = (sys.ugain + uncer_real);
    
    ftay = [                                                                                                                                                                                                                                              vktmax*cos(10*xorb3)*cos(xorb1)*cos(xorb2) - vktmax*cos(10*xorb3)*cos(xorb2)*sin(xorb1)*(x1 - xorb1) - vktmax*cos(10*xorb3)*cos(xorb1)*sin(xorb2)*(x2 - xorb2) - 10*vktmax*sin(10*xorb3)*cos(xorb1)*cos(xorb2)*(x3 - xorb3) - (vktmax*cos(10*xorb3)*cos(xorb1)*cos(xorb2)*(x1 - xorb1)^2)/2 - (vktmax*cos(10*xorb3)*cos(xorb1)*cos(xorb2)*(x2 - xorb2)^2)/2 - 50*vktmax*cos(10*xorb3)*cos(xorb1)*cos(xorb2)*(x3 - xorb3)^2 + (vktmax*cos(10*xorb3)*cos(xorb2)*sin(xorb1)*(x1 - xorb1)^3)/6 + (vktmax*cos(10*xorb3)*cos(xorb1)*sin(xorb2)*(x2 - xorb2)^3)/6 + (500*vktmax*sin(10*xorb3)*cos(xorb1)*cos(xorb2)*(x3 - xorb3)^3)/3 + vktmax*cos(10*xorb3)*sin(xorb1)*sin(xorb2)*(x1 - xorb1)*(x2 - xorb2) + 10*vktmax*sin(10*xorb3)*cos(xorb2)*sin(xorb1)*(x1 - xorb1)*(x3 - xorb3) + 10*vktmax*sin(10*xorb3)*cos(xorb1)*sin(xorb2)*(x2 - xorb2)*(x3 - xorb3) + (vktmax*cos(10*xorb3)*cos(xorb1)*sin(xorb2)*(x1 - xorb1)^2*(x2 - xorb2))/2 + (vktmax*cos(10*xorb3)*cos(xorb2)*sin(xorb1)*(x1 - xorb1)*(x2 - xorb2)^2)/2 + 50*vktmax*cos(10*xorb3)*cos(xorb2)*sin(xorb1)*(x1 - xorb1)*(x3 - xorb3)^2 + 5*vktmax*sin(10*xorb3)*cos(xorb1)*cos(xorb2)*(x1 - xorb1)^2*(x3 - xorb3) + 50*vktmax*cos(10*xorb3)*cos(xorb1)*sin(xorb2)*(x2 - xorb2)*(x3 - xorb3)^2 + 5*vktmax*sin(10*xorb3)*cos(xorb1)*cos(xorb2)*(x2 - xorb2)^2*(x3 - xorb3) - 10*vktmax*sin(10*xorb3)*sin(xorb1)*sin(xorb2)*(x1 - xorb1)*(x2 - xorb2)*(x3 - xorb3);...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               vktmax*sin(10*xorb3)*cos(xorb2) - (vktmax*sin(10*xorb3)*cos(xorb2)*(x2 - xorb2)^2)/2 - 50*vktmax*sin(10*xorb3)*cos(xorb2)*(x3 - xorb3)^2 + (vktmax*sin(10*xorb3)*sin(xorb2)*(x2 - xorb2)^3)/6 + 10*vktmax*cos(10*xorb3)*cos(xorb2)*(x3 - xorb3) - vktmax*sin(10*xorb3)*sin(xorb2)*(x2 - xorb2) - (500*vktmax*cos(10*xorb3)*cos(xorb2)*(x3 - xorb3)^3)/3 - 10*vktmax*cos(10*xorb3)*sin(xorb2)*(x2 - xorb2)*(x3 - xorb3) - 5*vktmax*cos(10*xorb3)*cos(xorb2)*(x2 - xorb2)^2*(x3 - xorb3) + 50*vktmax*sin(10*xorb3)*sin(xorb2)*(x2 - xorb2)*(x3 - xorb3)^2;...
            ugain*(x1 - xorb1)^3*((unom*cos(xorb2)*sin(xorb1))/6 - (cos(xorb1)*cos(xorb2)*(k1*Zi1_1 + k2*Zi2_1))/2) - ugain*(x2 - xorb2)*(unom*cos(xorb1)*sin(xorb2) - cos(xorb1)*cos(xorb2)*(k1*Zi1_2 + k2*Zi2_2)) - ugain*(x1 - xorb1)^2*((unom*cos(xorb1)*cos(xorb2))/2 + cos(xorb2)*sin(xorb1)*(k1*Zi1_1 + k2*Zi2_1)) - ugain*(x2 - xorb2)^2*((unom*cos(xorb1)*cos(xorb2))/2 + cos(xorb1)*sin(xorb2)*(k1*Zi1_2 + k2*Zi2_2)) - ugain*(x1 - xorb1)*(unom*cos(xorb2)*sin(xorb1) - cos(xorb1)*cos(xorb2)*(k1*Zi1_1 + k2*Zi2_1)) + ugain*(x2 - xorb2)^3*((unom*cos(xorb1)*sin(xorb2))/6 - (cos(xorb1)*cos(xorb2)*(k1*Zi1_2 + k2*Zi2_2))/2) + ugain*unom*cos(xorb1)*cos(xorb2) - ugain*(x1 - xorb1)*(x2 - xorb2)*(cos(xorb1)*sin(xorb2)*(k1*Zi1_1 + k2*Zi2_1) - unom*sin(xorb1)*sin(xorb2) + cos(xorb2)*sin(xorb1)*(k1*Zi1_2 + k2*Zi2_2)) + ugain*(x1 - xorb1)^2*(x2 - xorb2)*(sin(xorb1)*sin(xorb2)*(k1*Zi1_1 + k2*Zi2_1) + (unom*cos(xorb1)*sin(xorb2))/2 - (cos(xorb1)*cos(xorb2)*(k1*Zi1_2 + k2*Zi2_2))/2) + ugain*(x1 - xorb1)*(x2 - xorb2)^2*(sin(xorb1)*sin(xorb2)*(k1*Zi1_2 + k2*Zi2_2) + (unom*cos(xorb2)*sin(xorb1))/2 - (cos(xorb1)*cos(xorb2)*(k1*Zi1_1 + k2*Zi2_1))/2) + ugain*cos(xorb1)*cos(xorb2)*(x3 - xorb3)*(k1*Zi1_3 + k2*Zi2_3) - (ugain*cos(xorb1)*cos(xorb2)*(x1 - xorb1)^2*(x3 - xorb3)*(k1*Zi1_3 + k2*Zi2_3))/2 - (ugain*cos(xorb1)*cos(xorb2)*(x2 - xorb2)^2*(x3 - xorb3)*(k1*Zi1_3 + k2*Zi2_3))/2 - ugain*cos(xorb2)*sin(xorb1)*(x1 - xorb1)*(x3 - xorb3)*(k1*Zi1_3 + k2*Zi2_3) - ugain*cos(xorb1)*sin(xorb2)*(x2 - xorb2)*(x3 - xorb3)*(k1*Zi1_3 + k2*Zi2_3) + ugain*sin(xorb1)*sin(xorb2)*(x1 - xorb1)*(x2 - xorb2)*(x3 - xorb3)*(k1*Zi1_3 + k2*Zi2_3)];
end
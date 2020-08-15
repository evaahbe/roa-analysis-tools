%SCALED

function [ftay] = systemCL_taylor_o4_uncertain(xnom, unom, xloc, K,Pi,paras, uncertainty_vec)
    
    x1 = xloc(1);
    x2 = xloc(2);
    x3 = xloc(3);
    xnom1 = xnom(1);
    xnom2 = xnom(2);
    xnom3 = xnom(3);
    k1 = K(1);
    k2 = K(2);
    pi1_1 = Pi(1,1);
    pi1_2 = Pi(1,2);
    pi1_3 = Pi(1,3);
    pi2_1 = Pi(2,1);
    pi2_2 = Pi(2,2);
    pi2_3 = Pi(2,3);
    
    vktmax1 = paras.vw *paras.E0*1/paras.l + uncertainty_vec(1);
    vktmax2 = paras.vw *paras.E0*1/paras.l + uncertainty_vec(2);
    ugain = paras.vw*paras.E0*paras.Ktilde + uncertainty_vec(3);
    
    %see createTaylorForCopyPaste for derivation of taylor approx.
    %4th Order!!
    ftay = [                                                                                                                                                                                                                                              vktmax1*cos(10*xnom3)*cos(xnom1)*cos(xnom2) - vktmax1*cos(10*xnom3)*cos(xnom2)*sin(xnom1)*(x1 - xnom1) - vktmax1*cos(10*xnom3)*cos(xnom1)*sin(xnom2)*(x2 - xnom2) - 10*vktmax1*sin(10*xnom3)*cos(xnom1)*cos(xnom2)*(x3 - xnom3) - (vktmax1*cos(10*xnom3)*cos(xnom1)*cos(xnom2)*(x1 - xnom1)^2)/2 - (vktmax1*cos(10*xnom3)*cos(xnom1)*cos(xnom2)*(x2 - xnom2)^2)/2 - 50*vktmax1*cos(10*xnom3)*cos(xnom1)*cos(xnom2)*(x3 - xnom3)^2 + (vktmax1*cos(10*xnom3)*cos(xnom2)*sin(xnom1)*(x1 - xnom1)^3)/6 + (vktmax1*cos(10*xnom3)*cos(xnom1)*sin(xnom2)*(x2 - xnom2)^3)/6 + (500*vktmax1*sin(10*xnom3)*cos(xnom1)*cos(xnom2)*(x3 - xnom3)^3)/3 + vktmax1*cos(10*xnom3)*sin(xnom1)*sin(xnom2)*(x1 - xnom1)*(x2 - xnom2) + 10*vktmax1*sin(10*xnom3)*cos(xnom2)*sin(xnom1)*(x1 - xnom1)*(x3 - xnom3) + 10*vktmax1*sin(10*xnom3)*cos(xnom1)*sin(xnom2)*(x2 - xnom2)*(x3 - xnom3) + (vktmax1*cos(10*xnom3)*cos(xnom1)*sin(xnom2)*(x1 - xnom1)^2*(x2 - xnom2))/2 + (vktmax1*cos(10*xnom3)*cos(xnom2)*sin(xnom1)*(x1 - xnom1)*(x2 - xnom2)^2)/2 + 50*vktmax1*cos(10*xnom3)*cos(xnom2)*sin(xnom1)*(x1 - xnom1)*(x3 - xnom3)^2 + 5*vktmax1*sin(10*xnom3)*cos(xnom1)*cos(xnom2)*(x1 - xnom1)^2*(x3 - xnom3) + 50*vktmax1*cos(10*xnom3)*cos(xnom1)*sin(xnom2)*(x2 - xnom2)*(x3 - xnom3)^2 + 5*vktmax1*sin(10*xnom3)*cos(xnom1)*cos(xnom2)*(x2 - xnom2)^2*(x3 - xnom3) - 10*vktmax1*sin(10*xnom3)*sin(xnom1)*sin(xnom2)*(x1 - xnom1)*(x2 - xnom2)*(x3 - xnom3);...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               vktmax2*sin(10*xnom3)*cos(xnom2) - (vktmax2*sin(10*xnom3)*cos(xnom2)*(x2 - xnom2)^2)/2 - 50*vktmax2*sin(10*xnom3)*cos(xnom2)*(x3 - xnom3)^2 + (vktmax2*sin(10*xnom3)*sin(xnom2)*(x2 - xnom2)^3)/6 + 10*vktmax2*cos(10*xnom3)*cos(xnom2)*(x3 - xnom3) - vktmax2*sin(10*xnom3)*sin(xnom2)*(x2 - xnom2) - (500*vktmax2*cos(10*xnom3)*cos(xnom2)*(x3 - xnom3)^3)/3 - 10*vktmax2*cos(10*xnom3)*sin(xnom2)*(x2 - xnom2)*(x3 - xnom3) - 5*vktmax2*cos(10*xnom3)*cos(xnom2)*(x2 - xnom2)^2*(x3 - xnom3) + 50*vktmax2*sin(10*xnom3)*sin(xnom2)*(x2 - xnom2)*(x3 - xnom3)^2;...
            ugain*(x1 - xnom1)^3*((unom*cos(xnom2)*sin(xnom1))/6 - (cos(xnom1)*cos(xnom2)*(k1*pi1_1 + k2*pi2_1))/2) - ugain*(x2 - xnom2)*(unom*cos(xnom1)*sin(xnom2) - cos(xnom1)*cos(xnom2)*(k1*pi1_2 + k2*pi2_2)) - ugain*(x1 - xnom1)^2*((unom*cos(xnom1)*cos(xnom2))/2 + cos(xnom2)*sin(xnom1)*(k1*pi1_1 + k2*pi2_1)) - ugain*(x2 - xnom2)^2*((unom*cos(xnom1)*cos(xnom2))/2 + cos(xnom1)*sin(xnom2)*(k1*pi1_2 + k2*pi2_2)) - ugain*(x1 - xnom1)*(unom*cos(xnom2)*sin(xnom1) - cos(xnom1)*cos(xnom2)*(k1*pi1_1 + k2*pi2_1)) + ugain*(x2 - xnom2)^3*((unom*cos(xnom1)*sin(xnom2))/6 - (cos(xnom1)*cos(xnom2)*(k1*pi1_2 + k2*pi2_2))/2) + ugain*unom*cos(xnom1)*cos(xnom2) - ugain*(x1 - xnom1)*(x2 - xnom2)*(cos(xnom1)*sin(xnom2)*(k1*pi1_1 + k2*pi2_1) - unom*sin(xnom1)*sin(xnom2) + cos(xnom2)*sin(xnom1)*(k1*pi1_2 + k2*pi2_2)) + ugain*(x1 - xnom1)^2*(x2 - xnom2)*(sin(xnom1)*sin(xnom2)*(k1*pi1_1 + k2*pi2_1) + (unom*cos(xnom1)*sin(xnom2))/2 - (cos(xnom1)*cos(xnom2)*(k1*pi1_2 + k2*pi2_2))/2) + ugain*(x1 - xnom1)*(x2 - xnom2)^2*(sin(xnom1)*sin(xnom2)*(k1*pi1_2 + k2*pi2_2) + (unom*cos(xnom2)*sin(xnom1))/2 - (cos(xnom1)*cos(xnom2)*(k1*pi1_1 + k2*pi2_1))/2) + ugain*cos(xnom1)*cos(xnom2)*(x3 - xnom3)*(k1*pi1_3 + k2*pi2_3) - (ugain*cos(xnom1)*cos(xnom2)*(x1 - xnom1)^2*(x3 - xnom3)*(k1*pi1_3 + k2*pi2_3))/2 - (ugain*cos(xnom1)*cos(xnom2)*(x2 - xnom2)^2*(x3 - xnom3)*(k1*pi1_3 + k2*pi2_3))/2 - ugain*cos(xnom2)*sin(xnom1)*(x1 - xnom1)*(x3 - xnom3)*(k1*pi1_3 + k2*pi2_3) - ugain*cos(xnom1)*sin(xnom2)*(x2 - xnom2)*(x3 - xnom3)*(k1*pi1_3 + k2*pi2_3) + ugain*sin(xnom1)*sin(xnom2)*(x1 - xnom1)*(x2 - xnom2)*(x3 - xnom3)*(k1*pi1_3 + k2*pi2_3)];
end
function [sys,fns,pinit] = initializeKite(sys,fns,pinit,ana_type)


    load reference_UniSmoothScaled.mat
    
    sys.paras.l       = 60.00;         % line length
    sys.paras.vw      = 6.00;          % wind speed                  
    sys.paras.E0      = 5.68;          % glide ratio
    sys.paras.K       = 1;             % for scaled, otherwise 10.00; 
    sys.paras.vktmax  = sys.paras.vw  * sys.paras.E0 / sys.paras.l;
    sys.ugain         = sys.paras.vw * sys.paras.E0 * sys.paras.K / 8; %the constant with which u is multiplied to result in gammadot (times cos(theta)*cos(phi) missing to adjust for the wind direction)
    
    x_ref = reference_UniSmoothScaled{sys.paras.vw}(:,1:6);
    x = x_ref(:,1:3);
    u = x_ref(:,5);
    udot = x_ref(:,6);
    t = x_ref(:,4);
    tf = t(end);
    
    sys.xdim = 3;
    sys.udim = 1; 
    sys.Tperiod = tf;
    
    sys.xc = [0.3;0;0];
    
    sys.mu = 1;
    %sys.mu.nom = 1;
    %sys.mu.coeffs = [];

    pinit.Qmat = 0.1*eye(sys.xdim);
    pinit.Q0mat = 0.1*eye(sys.xdim);
    pinit.R = 1; 

    sys.rounds = 8; 

    fns.system_dyn= str2func('dynamics_Kite');
    fns.system_dyn_deriv = str2func('deriv_dynamics_Kite');
    fns.system_linearization = str2func('linear_dynamics_Kite');
    fns.system_dyn_poly = str2func('poly_dynamics_Kite');
    
    
    if strcmp(ana_type,'ROC')
        % for ROC computations provide the uncertain dynamics
        fns.system_dyn_poly_uncertain = str2func('uncertain_poly_dynamics_Kite');
        sys.uncermax = [0.4*sys.ugain]; %column vector of maximum variation of each uncertain parameter. If this is to be maximized just fill zeros or a desired minimum size
        sys.uncerdim = 1; % length of uncermax vector 
        [sys.permutvector,sys.uncer_count] = createUncertaintyPermutations(sys.uncerdim);
    end
    

    [sys.xorb, sys.uorb, sys.uorbdot, sys.tvec] = construct_spline_control(sys,fns,t,x,u,udot);
    [sys.eta, sys.v, sys.vdot] = hyperplane_base_control(sys,fns); %returns splines
    
end
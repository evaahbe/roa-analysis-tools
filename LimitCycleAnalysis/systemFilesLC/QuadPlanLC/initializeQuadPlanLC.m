function [sys,fns,pinit] = initializeQuadPlanLC(sys,fns,pinit,ana_type)


    load LC_QuadPlanLC_mu8.mat
    
    sys.xdim = 2;
    sys.udim = 0; 
    
    t = x_ref(:,1);
    x = x_ref(:,2:sys.xdim+1);
    %u = x_ref(:,sys.xdim+2:sys.xdim+2+sys.udim-1);
    %u = x_ref(:,sys.xdim+2+sys.udim:sys.xdim+2+sys.udim+sys.udim-1);
    
    tf = t(end);   
    
    sys.Tperiod = tf;
    
    sys.xc = [0;0];
    
    sys.mu = 8; %8 is nominal for quadplanlc

    pinit.Qmat = 0.1*eye(sys.xdim);
    pinit.Q0mat = 0.1*eye(sys.xdim);

    sys.rounds = 8; 

    fns.system_dyn= str2func('dynamics_quadplan');
    fns.system_dyn_deriv = str2func('deriv_dynamics_quadplan');
    fns.system_linearization = str2func('linear_dynamics_quadplan');
    fns.system_dyn_poly = str2func('dynamics_quadplan'); % dynamics for VDP already polynomial
    
    if strcmp(ana_type,'ROC')
        % for ROC computations provide the uncertain dynamics
        fns.system_dyn_poly_uncertain = str2func('uncertain_dynamics_quadplan');
        sys.uncermax = [0.01]; %column vector of maximum variation of each uncertain parameter. If this is to be maximized just fill zeros or a desired minimum size
        sys.uncerdim = 1; % length of uncermax vector, use 0 if no uncertainty! 
        [sys.permutvector,sys.uncer_count] = createUncertaintyPermutations(sys.uncerdim);
    end

    if sys.udim ==0
        [sys.xorb, sys.tvec] = construct_spline_auto(sys,fns,t,x);
        [sys.eta, sys.v, sys.vdot] = hyperplane_base_auto(sys,fns); %returns splines
    else
        pinit.R = 1; 
        [sys.xorb,sys.uorb,sys.uorbdot, sys.tvec] = construct_spline_control(sys,fns,t,x,u,udot);
        [sys.eta, sys.v, sys.vdot] = hyperplane_base_control(sys,fns); %returns splines
    end
end
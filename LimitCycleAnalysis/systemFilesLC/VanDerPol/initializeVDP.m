function [sys,fns,pinit] = initializeVDP(sys,fns,pinit,ana_type)


    load VDP_LCtraj.mat
    % use VDP_LCtraj.mat for less accuracy but significant speedup (for
    % Lyapunov equation )
    
    x = xtraj{2};
    t = xtraj{1};
    
    tf = t(end);
    
    sys.xdim = 2;
    sys.udim = 0; 
    sys.Tperiod = tf;
    
    sys.xc = [0.;0.];
    
    sys.mu = 1;
    %sys.mu.nom = 1;
    %sys.mu.coeffs = [];

    pinit.Qmat = 0.1*eye(sys.xdim);
    pinit.Q0mat = 0.1*eye(sys.xdim);

    sys.rounds = 8; 

    fns.system_dyn= str2func('dynamics_vdp');
    fns.system_dyn_deriv = str2func('deriv_dynamics_vdp');
    fns.system_linearization = str2func('linear_dynamics_vdp');
    fns.system_dyn_poly = str2func('dynamics_vdp'); % dynamics for VDP already polynomial
    
    if strcmp(ana_type,'ROC')
        % for ROC computations provide the uncertain dynamics
        fns.system_dyn_poly_uncertain = str2func('uncertain_dynamics_vdp');
        sys.uncermax = [0.2*sys.mu]; %column vector of maximum variation of each uncertain parameter. If this is to be maximized just fill zeros or a desired minimum size
        sys.uncerdim = 1; % length of uncermax vector 
        [sys.permutvector,sys.uncer_count] = createUncertaintyPermutations(sys.uncerdim);
    end

    [sys.xorb, sys.tvec] = construct_spline_auto(sys,fns,t,x);
    [sys.eta, sys.v, sys.vdot] = hyperplane_base_auto(sys,fns); %returns splines
    
end
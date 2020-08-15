%%% THIS IS THE TEMPLATE FILE FOR A BOTH ROA AND ROC ANALYSIS OF A LIMIT
%%% CYCLE 

function [sys,fns,pinit] = initializeTemplate(sys,fns,pinit,ana_type)

    % NOTE: if any changes are made to the dynamics (e.g. parameter change)
    % the Lyapunov initial solution file has to be manually deleted in the 
    % folder 'initialLyapunov' 

    %%%====================================================================
    %%%====================================
    % CHANGE ALL OF THE BELOW ACCORDINGLY:
    
    load TEMPLATE_LC_trajectory.mat
    
    sys.xdim = 2;   % number of states
    sys.udim = 0;   % number of inputs
    
    t = x_ref(:,1);             % trajectory time vector. This needs to be according to the LC trajectory data provided in the loaded file above. Change if necessary.
    x = x_ref(:,2:sys.xdim+1);  % trajectory state matrix. This needs to be according to the LC trajectory data provided in the loaded file above. Change if necessary.
    %u = x_ref(:,sys.xdim+2:sys.xdim+2+sys.udim-1);                         % trajectory input matrix (if applicable)
    %udot = x_ref(:,sys.xdim+2+sys.udim:sys.xdim+2+sys.udim+sys.udim-1);    % trajectory input derivative matrix (if available)
    
    sys.xc = [0;0];                     % center point location if cp-MOC chosen (sys.trafo = 2 in main file)
    
    sys.mu(1) = 1;                      % parameter vector in case desired, for easy testing of different dynamics

    pinit.Qmat  = 0.1*eye(sys.xdim);    % pos-def. matrix for Lyapunov/Riccati equation
    pinit.Q0mat = 0.1*eye(sys.xdim);    % pos-def. initial matrix for Lyapunov/Riccati equation

    sys.rounds = 8;                     % number of periods for convergence of Lyapunov equation solution (last period is used)

    %replace 'TEMPLATE' with actual function name
    fns.system_dyn= str2func('dynamics_TEMPLATE');                      % true system dynamics
    fns.system_dyn_deriv = str2func('deriv_dynamics_TEMPLATE');         % dynamic derivative of the true dynamics
    fns.system_linearization = str2func('linear_dynamics_TEMPLATE');    % Jacobian file
    fns.system_dyn_poly = str2func('poly_dynamics_TEMPLATE');           % polynomial approx. of dynamics (or true ones, if already poly)              
    
    
    % only specify if uncertain dynamics are considered (ROC analysis)
    sys.uncermax = [0.*sys.mu];     % column vector of maximum absolute variation of each uncertain parameter. If this is to be maximized just fill very small values (e.g., 1e-3) or desired values
    fns.system_dyn_poly_uncertain = str2func('uncertain_dynamics_TEMPLATE'); % uncertain dynamics in polynomial approx. (if not poly already)
    
    % UNTIL HERE
    %%%====================================
    %%%====================================================================
    
    
    sys.Tperiod = t(end);
    
    if strcmp(ana_type,'ROC')
        sys.uncerdim = sum(sys.uncermax~=0);               
        [sys.permutvector,sys.uncer_count] = createUncertaintyPermutations(sys.uncerdim);
    end

    if sys.udim ==0
        [sys.xorb, sys.tvec] = construct_spline_auto(sys,fns,t,x);
        [sys.eta, sys.v, sys.vdot] = hyperplane_base_auto(sys,fns);
    else
        pinit.R = 1; 
        [sys.xorb,sys.uorb,sys.uorbdot, sys.tvec] = construct_spline_control(sys,fns,t,x,u,udot);
        [sys.eta, sys.v, sys.vdot] = hyperplane_base_control(sys,fns); 
    end
end
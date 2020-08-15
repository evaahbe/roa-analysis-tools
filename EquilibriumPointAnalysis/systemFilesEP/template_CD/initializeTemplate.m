%%% THIS IS THE TEMPLATE FILE FOR A FEEDBACK CONTROLLED STOCHASTIC SYSTEM
%%% WITH AN EQUILIBRIUM POINT

function [sys,fns] = initializeTemplate(sys,fns)
    
    %%%====================================================================
    %%%====================================
    % CHANGE ALL OF THE BELOW ACCORDINGLY:
    
    % Replace TEMPLATE with system file name (also in the function name above!)
    fns.system_dyn= str2func('dynamics_TEMPLATE');            % true system dynamics
    fns.system_dyn_int= str2func('dynamics_TEMPLATE_INT');    % true dynamics with time para, needed for EP computation
    fns.system_dyn_poly = str2func('poly_dynamics_TEMPLATE'); % polynomial approximation (or true system, if polynomial)
    
    sys.xdim = 2; % number of states
    sys.udim = 1; % number of inputs
    
    % Uncertainty type (has to be either or)
     sys.distType = 'uniform'; %'normal'
    
    % uncomment and enter according to distribution type 
    %  - for uniform distribution:
    sys.mu{1}.mu_A = 0; % lower distribution limit of first random variable
    sys.mu{1}.mu_B = 1; % upper distribution limit of first random variable
     
    % - for normal distribution:
%     sys.mu{1}.mu_A = 1; % lower distribution limit
%     sys.mu{1}.mu_B = 0.1; % upper distribution limit
    
    sys.varfix = 0*eye(sys.xdim); %fixed covariance matrix for initial states in stochastic ROA

    hx = sym('hx',[1 sys.xdim]); %Don't change this.
    
    % enter desired state feedback vector:
    sys.h_vec = [hx(1);hx(2)^2]; %hx(1) = first state, hx(2) = second state,...
    
    % Optional: Input constraints (IC)
    sys.include_input_con = 1; %-2 = nominal, -1 = stochastic OL, 0 = CL without IC, 1 = CL with IC
    
    sys.umin = -0.5; %lower input limit (only relevant if sys.include_input_con = 1)
    sys.umax = 0.5; %upper input limit (only relevant if sys.include_input_con = 1)
    
    % UNTIL HERE
    %%%====================================
    %%%====================================================================
    
    if strcmp(sys.distType,'uniform')
        for i = 1:length(sys.mu)
            sys.mu{1}.mu_nom = sys.mu{1}.mu_A; 
            sys.mu{1}.mu_coefs = galerkin_projection_normal(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p); 
        end
    elseif strcmp(sys.distType, 'normal')
        for i = 1:length(sys.mu)
            sys.mu{2}.mu_nom = (sys.mu{2}.mu_A+sys.mu{2}.mu_B)/2; 
            sys.mu{2}.mu_coefs = galerkin_projection_uniform(sys.mu{2}.mu_A,sys.mu{2}.mu_B,sys.p); 
        end
    end
     
    if sys.include_input_con == -2 %sets all to zero for nominal case
        sys.p = 0;
        sys.mu{1}.mu_A = sys.mu{1}.mu_nom;
        sys.mu{1}.mu_B = sys.mu{1}.mu_nom;
        sys.mu{1}.mu_coefs = galerkin_projection_uniform(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p);
        sys.mu{2}.mu_A = sys.mu{2}.mu_nom;
        sys.mu{2}.mu_B = sys.mu{2}.mu_nom;
        sys.mu{2}.mu_coefs = galerkin_projection_uniform(sys.mu{2}.mu_A,sys.mu{2}.mu_B,sys.p);
    end
    
    
end

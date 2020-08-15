%%% THIS IS THE TEMPLATE FILE FOR A STOCHASTIC SYSTEM
%%% WITH AN EQUILIBRIUM POINT


function [sys,fns] = initializeTemplate(sys,fns)
    
    %%%====================================================================
    %%%====================================
    % CHANGE ALL OF THE BELOW ACCORDINGLY:
    
    % Replace TEMPLATE with system file name (also in the function name above!)
    fns.system_dyn= str2func('dynamics_TEMPLATE');              % true system dynamics
    fns.system_dyn_int= str2func('dynamics_TEMPLATE_INT');      % true dynamics with time para, needed for EP computation
    fns.system_dyn_poly = str2func('poly_dynamics_TEMPLATE');   % polynomial approximation (or true system, if polynomial)
    
    sys.xdim = 2; % number of states
    sys.udim = 0; % number of inputs
     
     % Uncertainty type (has to be either or)
     sys.distType = 'uniform'; % other option: 'normal'
    
    % uncomment and enter according to distribution type 
    %  - for uniform distribution:
    sys.mu{1}.mu_A = 0; % lower distribution limit of first random variable
    sys.mu{1}.mu_B = 1; % upper distribution limit of first random variable
     
    % - for normal distribution:
%     sys.mu{1}.mu_A = 1;   % lower distribution limit
%     sys.mu{1}.mu_B = 0.1; % upper distribution limit
    
    sys.varfix = 0*eye(sys.xdim); %fixed covariance matrix for initial states in stochastic ROA
    
    % UNTIL HERE
    
    % Change these only if needed
    sys.Tintegral = 100;         % integration time for PCE EP simulation (function findPCEEPsimu) 
    sys.x0 = zeros(sys.xdim,1);  % initial condition for PCE EP simulation (function findPCEEPsimu)
    
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
  
end

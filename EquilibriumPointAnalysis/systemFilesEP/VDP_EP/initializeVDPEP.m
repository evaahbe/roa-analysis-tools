function [sys,fns] = initializeVDPEP(sys,fns)
    
    sys.xdim = 2;
    sys.udim = 0; 
    
    sys.Tintegral = 50;         % integration time for PCE EP simulation (function findPCEEPsimu) 
    sys.x0 = zeros(sys.xdim,1); % initial condition for PCE EP simulation (function findPCEEPsimu)
    
    sys.distType = 'uniform'; %'normal'
    
    % enter according to distribution type
     
    %for uniform distribution
    sys.mu{1}.mu_A = 0.7; % lower distribution limit
    sys.mu{1}.mu_B = 1.3; % upper distribution bound
    sys.mu{1}.mu_nom = (sys.mu{1}.mu_A+sys.mu{1}.mu_B)/2; % nominal mu value;
    sys.mu{1}.mu_coefs = galerkin_projection_uniform(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p); 
    
    % for normal distribution
%     sys.mu{1}.mu_A = 1; % lower distribution limit
%     sys.mu{1}.mu_B = 0.1; % upper distribution bound
%     sys.mu{1}.mu_nom = sys.mu{1}.mu_A+; % nominal mu value;
%     sys.mu{1}.mu_coefs = galerkin_projection_normal(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p); 


    %sys.varfix = 0*eye(sys.xdim);
    sys.varfix = 0.0*eye(sys.xdim);
    %sys.varfix{3} = 0.03*eye(sys.xdim);
    
    fns.system_dyn= str2func('dynamics_vdp_ep');
    fns.system_dyn_int= str2func('dynamics_vdp_ep_INT');
    fns.system_dyn_poly = str2func('dynamics_vdp_ep'); % dynamics for VDP already polynomial

    
end
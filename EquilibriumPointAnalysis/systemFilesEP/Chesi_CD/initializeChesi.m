function [sys,fns] = initializeChesi(sys,fns)
    
    sys.xdim = 2;
    sys.udim = 1; 
    
    sys.Tintegral = 100;         % integration time for PCE EP simulation (function findPCEEPsimu) 
    sys.x0 = zeros(sys.xdim,1); % initial condition for PCE EP simulation (function findPCEEPsimu)
    
    sys.distType = 'uniform'; %'normal'
    
    % enter according to distribution type
     
    %for uniform distribution
    sys.mu{1}.mu_A = 0.8; % lower distribution limit
    sys.mu{1}.mu_B =1.2; % upper distribution bound
    sys.mu{1}.mu_nom = (sys.mu{1}.mu_A+sys.mu{1}.mu_B)/2; % nominal mu value;
    sys.mu{1}.mu_coefs = galerkin_projection_uniform(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p); 
%     sys.mu{2}.mu_A = 0.5; % lower distribution limit
%     sys.mu{2}.mu_B = 1.5; % upper distribution bound
%     sys.mu{2}.mu_nom = (sys.mu{2}.mu_A+sys.mu{2}.mu_B)/2; % nominal mu value;
%     sys.mu{2}.mu_coefs = galerkin_projection_uniform(sys.mu{2}.mu_A,sys.mu{2}.mu_B,sys.p); 
    
    sys.varfix = 0*eye(sys.xdim);
  
    % for normal distribution
%     sys.mu{1}.mu_A = 1; % lower distribution limit
%     sys.mu{1}.mu_B = 0.1; % upper distribution bound
%     sys.mu{1}.mu_nom = sys.mu{1}.mu_A+; % nominal mu value;
%     sys.mu{1}.mu_coefs = galerkin_projection_normal(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p); 

    hx = sym('hx',[1 sys.xdim]); 
    % enter desired state feedback vector:
    sys.h_vec = [hx(1);hx(2)]; %hx(1) = first state, hx(2) = second state,...
    %sys.h_vec = [hx(1);hx(2);hx(1)^2;hx(1)*hx(2);hx(2)^2];
    
    % Optional: Input constraints
    sys.include_input_con = 0; %-2 = nominal, -1 = stochastic OL, 0 = CL without IC, 1 = CL with IC
    sys.umin = -0.5;
    sys.umax = 0.5;
    
    if sys.include_input_con == -2
        sys.p = 0;
        sys.mu{1}.mu_A = sys.mu{1}.mu_nom;
        sys.mu{1}.mu_B = sys.mu{1}.mu_nom;
        sys.mu{1}.mu_coefs = galerkin_projection_uniform(sys.mu{1}.mu_A,sys.mu{1}.mu_B,sys.p);
    end
    

    
    fns.system_dyn= str2func('dynamics_chesi');
    fns.system_dyn_int= str2func('dynamics_chesi_INT');
    fns.system_dyn_poly = str2func('dynamics_chesi'); % dynamics for VDP already polynomial

    
end

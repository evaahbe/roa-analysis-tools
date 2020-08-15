function veriROA_PCE(sys,fns,filenames,numsets)

    nr_mu = length(sys.mu);
    mu = sdpvar(nr_mu,1);
    
    x = sdpvar(sys.xdim,1);
    xi = sdpvar((sys.p+1)*sys.xdim,1);
    
    % PCE preprocessing
    vartot = [mu;x];
    
    varitot = [];
    for j = 1:nr_mu
        varitot = [varitot;sys.mu{j}.mu_coefs];
    end
    for j = 1:sys.xdim
        varitot = [varitot;xi((sys.p+1)*(j-1)+1:(sys.p+1)*j)'];
    end
    
    fxstoch = fns.system_dyn_poly(x,mu);
    fxi_preshift = castPCEdynamics(vartot,varitot,sys,fxstoch);

    % get PCE equilibrium point via simulation of PCE system
    try             % load results of already performed iteration steps
        load(filenames.intermedresultsROAPCE)
    catch 
        xiep = findPCEEPsimu(xi,fxi_preshift,sys)';
        Qvals = [];
        save(filenames.intermedresultsROAPCE,'xiep','Qvals')
    end
    
    %coordinate shift:
    fxi = replace(fxi_preshift,xi,xi+xiep);
    inicalc.fxi = clean(fxi.',numsets.clean_thresh);

    if isempty(Qvals) 
        % linearization around zero point and Lyapunov equation result
        Axiep = replace(clean(jacobian(fxi.',xi),numsets.clean_thresh),xi,zeros(size(xi))); 
        Qmat = diag([ones(size(Axiep,1),1)]);
        Pmat = lyap(transpose(Axiep),Qmat)/numsets.initVscale;
        inicalc.V = xi.'*Pmat*xi;
        inicalc.dVdx = clean(jacobian(inicalc.V,xi),numsets.clean_thresh);
        multis = [];
        Smat = 1.5*Pmat;
        Smatdetlist = [10*det(Smat),det(Smat)];
        gamma   = 1;  % V sublevel set,  fixed
        alpha   = 1;  % B sublevel set,  fixed
        save(filenames.intermedresultsROAPCE,'-append','Smatdetlist','gamma','alpha','multis','Smat','xiep')
    end


    gterm = numsets.gfac*xi'*xi;
    degs = numsets.degs;

    for jj=length(Smatdetlist)-1:numsets.iteration_max
        
        if Smatdetlist(end-1)- Smatdetlist(end) < numsets.convCrit*Smatdetlist(end)
            break
        end
        fprintf(1,'\n \n This is iteration nr.: %d \n',jj)
        
        multis = find_multiplier(xi,numsets,degs,Qvals,inicalc,gterm,Smat,Smatdetlist,multis);
        if multis.flag ==-1
            break
        else
            save(filenames.intermedresultsROAPCE,'-append','multis')
        end
        [Smat,Qvals] = optimize_V(xi,numsets,degs,inicalc,multis,Smat,gterm,Smatdetlist,Qvals);
        if Smat==-1
            Smat = Smatdetlist(end);
            break
        else
            Smatdetlist = [Smatdetlist,det(Smat)];
            save(filenames.intermedresultsROAPCE,'-append','Smatdetlist','Qvals','Smat')
        end
    end
    save(filenames.finalresultsROAPCE,'Qvals','multis','Smat','Smatdetlist','gamma','alpha','xiep','sys','numsets')
end


function multis = find_multiplier(xi,numsets,degs,Qvals,inicalc,gterm,Smat,Smatlist,multisold)

    alpha = 1;
    gamma = 1;
                 
    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);
    
    [s1,c1,~] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [s2,c2,~] = polynomial([xi],degs.s2_dU,degs.s2_dL);
    
 
    if isempty(Qvals)           
        V = inicalc.V; 
        dVdx = inicalc.dVdx;
    else   
        V = mV.'*Qvals*mV; 
        dVdx = clean(jacobian(V,xi), numsets.clean_thresh);        
    end
    
    dVdt = dVdx*inicalc.fxi;
    
    eq1 = - dVdt - s1 * (alpha-V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2 + (V - alpha));
           
    F1 = [sos(s1),sos(eq1)];
    F2 = [sos(s2),sos(eq2)];
    
    [sol1] = solvesos(F1,[],numsets.sdpsetting1,c1); 
    [sol2] = solvesos(F2,[],numsets.sdpsetting1,c2);

    res1 = sol1.problem;
    res2 = sol2.problem;
    
    if (res1==-1 || res2==-1) && length(Smatlist)>=5
        multis.flag =-1;
        multis.c1 = multisold.c1; 
        multis.c2 = multisold.c2;
        fprintf(1,'-------Ran into numerical problems, final results file created but terminating early..-------\n\n')
        return
    elseif res1~=0 || res2~=0
        fprintf(1,'The following errors were found: %d %d \n',[res1,res2]);
        error('Error in multiplier search')
    end
       
    multis.c1 = value(c1)'; 
    multis.c2 = value(c2)';
    multis.flag = 0;
                                   
    fprintf(1,'-------Multiplier search successfully finished-------\n\n')

end


function [Smatres,Qvalres] = optimize_V(xi,numsets,degs,inicalc,multis,Smatold,gterm,Smatlist,Qvalsold)

    alpha = 1;
    gamma = 1;
    
    Smatdiff = sdpvar(length(xi),length(xi));
    Smat = Smatold - Smatdiff;
    
    [~,~,m1] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([xi],degs.s2_dU,degs.s2_dL);
   
    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);
       
    Qvals = sdpvar(size(mV,1),size(mV,1)); 
    V = mV.'*Qvals*mV;      
    dVdx = clean(jacobian(V,xi), numsets.clean_thresh);
    dVdt = dVdx*inicalc.fxi;
               
    s1fix = multis.c1*m1;
    s2fix = multis.c2*m2;
    
    eq1 = - dVdt - s1fix * (alpha - V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2fix + (V - alpha)); 

    %F = [sos(eq1), sos(eq2), sos(mV.'*Qvals*mV - gterm), sos(xi.'*Smat*xi)]; 
    F = [sos(eq1), sos(eq2), Qvals>=0, Smat>=0,Smatdiff>=0];
    
    [sol] = solvesos(F, -geomean(Smatdiff),numsets.sdpsetting2,[Qvals(:);Smatdiff(:)]);
    
    if sol.problem == -1 && length(Smatlist)>=5
        Smatres = -1;
        Qvalres = Qvalsold;
        fprintf(1,'-------Ran into numerical problems, final results file created but terminating early..-------\n\n')
        return
    elseif sol.problem ~= 0 
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in V optimization')
    end
    
    Smatres = value(Smat);
    Qvalres = value(Qvals);
    fprintf('This is the geomean result: %f \n',det(Smatres));  
     
    fprintf(1,'-------V optimization successfully finished-------\n\n')
end


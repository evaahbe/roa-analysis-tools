function veriROA_PCE_CD(sys,fns,filenames,numsets)

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
    
    h_vec = make_hvec_sdpvar(x,sys);
    
    K = sdpvar(sys.udim,length(h_vec));
    u = K*h_vec;
    Kvals = zeros(size(K));
    if length(numsets.initKvals(:)) ==1 %scalar initial value
        initKvals = numsets.initKvals * ones(size(K));
    else                                %initial value with right size
        initKvals = numsets.initKvals .* ones(size(K));
    end
    
    if length(numsets.initt(:)) ==1
        initt = numsets.initt*ones(sys.udim,sys.p);
    else
        initt = numsets.initt.*ones(sys.udim,sys.p);
    end
  
    tvar = [];
        
    fxstoch = fns.system_dyn_poly(x,mu,u);
    fxi_preshift = castPCEdynamics(vartot,varitot,sys,fxstoch);
    
    fxi_for_ep = replace(fxi_preshift,K,Kvals);

    % get PCE equilibrium point via simulation of PCE system
    try             % load results of already performed iteration steps
        load(filenames.intermedresultsROAPCE_pre,'xiep')
    catch 
        xiep = findPCEEPsimu(xi,fxi_for_ep,sys)';
    end
    
    %coordinate shift:
    fxi = replace(fxi_preshift,xi,xi+xiep);
    inicalc.fxi = clean(fxi.',numsets.clean_thresh);
    
    normsvec(1) = 0;
    if strcmp(sys.distType,'uniform')
        for i = 2:sys.p+1
            normsvec(i) = legendreProduct(i-1,i-1);
        end
    elseif strcmp(sys.distType,'normal')
        for i = 2:sys.p+1
            normsvec(i) = hermiteProduct(i-1,i-1);
        end
    end
    
    ulawall = castPCEdynamics(vartot, varitot,sys,u);
    for i = 1:sys.udim
       ulaw.u0(i) = ulawall((i-1)*(sys.p+1)+1);
       ulaw.ucomp(i,:) = sqrt(normsvec(2:end)).*ulawall((i-1)*(sys.p+1)+2:i*(sys.p+1));
    end
    
    init_flag = 1;
    veri_step = 1; % 1 for the pre ROA calc, 2 for ROA with control design
    
    if exist(filenames.intermedresultsROAPCE,'file')
        load(filenames.intermedresultsROAPCE)   
        if length(Smatdetlist_fb)<=2
            Kvals = initKvals; %in this case the init value is possibily being tuned, allow for the new one to be assigned
        end
        veri_step = 2;
    elseif exist(filenames.finalresultsROAPCE_pre,'file')
        load(filenames.finalresultsROAPCE_pre,'Smatdetlist','Qvals','gamma','alpha','multis','Smat')  
        veri_step = 2;
        Smatdetlist_fb = [10*Smatdetlist(end),Smatdetlist(end)];
        Kvals = initKvals; % initialize Kvals for the FB design loop with very small values (to be tuned)
        save(filenames.intermedresultsROAPCE,'Smatdetlist_fb','Qvals','gamma','alpha','multis','Smat','xiep','Kvals')
    elseif exist(filenames.intermedresultsROAPCE_pre,'file')
        load(filenames.intermedresultsROAPCE_pre,'Smatdetlist','Qvals','gamma','alpha','multis','Smat')
        if isempty(Qvals) 
            init_flag = 0;
        end
    else
        init_flag = 0;
    end
    if init_flag == 0
        % linearization around zero point and Lyapunov equation result
        Axiep = replace(clean(jacobian(fxi.',xi),numsets.clean_thresh),[xi;K(:)],[zeros(length(xi),1);Kvals(:)]); 
        Qmat = diag([ones(size(Axiep,1),1)]);
        Pmat = lyap(transpose(Axiep),Qmat)/numsets.initVscale;
        inicalc.V = xi.'*Pmat*xi;
        inicalc.dVdx = clean(jacobian(inicalc.V,xi),numsets.clean_thresh);
        Qvals = [];
        multis = [];     
        Smat = 1.5*Pmat;
        Smatdetlist = [1e30,det(Smat)];
        gamma   = 1;  % V sublevel set, fixed
        alpha   = 1;  % B sublevel set, fixed
        save(filenames.intermedresultsROAPCE_pre,'Smatdetlist','Qvals','gamma','alpha','multis','Smat','xiep')
    end

    gterm = numsets.gfac*xi'*xi;
    degs = numsets.degs;
   
    if veri_step == 1
        
        fprintf(1, 'Running the pre ROA calculations first...')
        for jj=(length(Smatdetlist)-1):numsets.iteration_max_preROA
        
            if Smatdetlist(end-1)- Smatdetlist(end) < numsets.convCrit_preROA*Smatdetlist(end)
                break
            end
            fprintf(1,'\n \n This is pre ROA iteration nr.: %d \n',jj)
        
            multis = find_multiplier(xi,K,numsets,degs,Qvals,inicalc,gterm,Smat,Smatdetlist,multis,Kvals);
            if multis.flag ==-1
                break
            else
                save(filenames.intermedresultsROAPCE_pre,'-append','multis')
            end
            [Smat,Qvals] = optimize_V(xi,K,numsets,degs,inicalc,multis,Smat,gterm,Smatdetlist,Qvals,Kvals);
            if Smat==-1
                break
            else
                Smatdetlist = [Smatdetlist,det(Smat)];
                save(filenames.intermedresultsROAPCE_pre,'-append','Smatdetlist','Qvals','Smat')
            end
        end
        save(filenames.finalresultsROAPCE_pre,'Qvals','multis','Smat','Smatdetlist','gamma','alpha','xiep','sys','numsets')
        Smatdetlist_fb = [10*Smatdetlist(end),Smatdetlist(end)];
        Kvals = initKvals;
        save(filenames.intermedresultsROAPCE,'Smatdetlist_fb','Qvals','gamma','alpha','multis','Smat','xiep','Kvals')
        fprintf(1, 'Pre ROA calculations completed.')
    end
    
    if sys.include_input_con>=-0.01 % -1 for stochastic OL and -2 for nominal OL
    
        fprintf(1, 'Running the control design for maximal ROA calculations...')

        if sys.include_input_con ==1

            if isempty(tvar)
                tvar = initt;
            end

            for jj = (length(Smatdetlist_fb)-1):numsets.iteration_max

                if Smatdetlist_fb(end-1)- Smatdetlist_fb(end) < numsets.convCrit*Smatdetlist_fb(end)
                    break
                end
                fprintf(1,'\n \n This is FB ROA iteration nr.: %d \n',jj)

                multis = find_multiplier_wIC(xi,K,sys,numsets,degs,Qvals,inicalc,gterm,Smat,Smatdetlist_fb,multis,Kvals,tvar,ulaw);
                if multis.flag ==-1
                    break
                else
                    save(filenames.intermedresultsROAPCE,'-append','multis')
                end

                [Kvals,tvar] = optimize_K_wIC(xi,K,sys,numsets,degs,inicalc,multis,Smat,gterm,Qvals,ulaw);
                save(filenames.intermedresultsROAPCE,'-append','Kvals','tvar')

                [Smat,Qvals] = optimize_V_wIC(xi,K,sys,numsets,degs,inicalc,multis,Smat,gterm,Smatdetlist_fb,Qvals,Kvals,tvar,ulaw);

                if Smat==-1
                    Smat = Smatdetlist_fb(end);
                    break
                else
                    Smatdetlist_fb = [Smatdetlist_fb,det(Smat)];
                    save(filenames.intermedresultsROAPCE,'-append','Smatdetlist_fb','Qvals','Smat')
                end
            end
            save(filenames.finalresultsROAPCE,'Qvals','multis','Smat','Smatdetlist_fb','gamma','alpha','xiep','sys','numsets','Kvals','tvar')
        else
            for jj = (length(Smatdetlist_fb)-1):numsets.iteration_max

                if Smatdetlist_fb(end-1)- Smatdetlist_fb(end) < numsets.convCrit*Smatdetlist_fb(end)
                    break
                end
                fprintf(1,'\n \n This is FB ROA iteration nr.: %d \n',jj)

                multis = find_multiplier(xi,K,numsets,degs,Qvals,inicalc,gterm,Smat,Smatdetlist_fb,multis,Kvals);
                if multis.flag ==-1
                    break
                else
                    save(filenames.intermedresultsROAPCE,'-append','multis')
                end

                [Kvals] = optimize_K(xi,K,numsets,degs,Qvals,inicalc,gterm,multis);
                save(filenames.intermedresultsROAPCE,'-append','Kvals')

                [Smat,Qvals] = optimize_V(xi,K,numsets,degs,inicalc,multis,Smat,gterm,Smatdetlist_fb,Qvals,Kvals);

                if Smat==-1
                    Smat = Smatdetlist_fb(end);
                    break
                else
                    Smatdetlist_fb = [Smatdetlist_fb,det(Smat)];
                    save(filenames.intermedresultsROAPCE,'-append','Smatdetlist_fb','Qvals','Smat')
                end
            end
            save(filenames.finalresultsROAPCE,'Qvals','multis','Smat','Smatdetlist_fb','gamma','alpha','xiep','sys','numsets','Kvals')         
        end
    else
        
        for jj=(length(Smatdetlist)-1):numsets.iteration_max
        
            if Smatdetlist(end-1)- Smatdetlist(end) < numsets.convCrit*Smatdetlist(end)
                break
            end
            fprintf(1,'\n \n This is pre ROA iteration nr.: %d \n',jj)
        
            multis = find_multiplier(xi,K,numsets,degs,Qvals,inicalc,gterm,Smat,Smatdetlist,multis,Kvals);
            if multis.flag ==-1
                break
            else
                save(filenames.intermedresultsROAPCE_pre,'-append','multis')
            end
            [Smat,Qvals] = optimize_V(xi,K,numsets,degs,inicalc,multis,Smat,gterm,Smatdetlist,Qvals,Kvals);
            if Smat==-1
                break
            else
                Smatdetlist = [Smatdetlist,det(Smat)];
                save(filenames.intermedresultsROAPCE_pre,'-append','Smatdetlist','Qvals','Smat')
            end
        end
        save(filenames.finalresultsROAPCE,'Qvals','multis','Smat','Smatdetlist','gamma','alpha','xiep','sys','numsets')
    end
        
       
      
    
end


function multis = find_multiplier_wIC(xi,K,sys,numsets,degs,Qvals,inicalc,gterm,Smat,Smatlist,multisold,Kvals,tvar,ulaw)
    
    gamma = 1;
    alpha = 1;
  
    mV = monolist([xi],degs.V_dU/2,degs.V_dL/2);

    [s1,c1,~] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [s2,c2,~] = polynomial([xi],degs.s2_dU,degs.s2_dL);

    
    fxi = replace(inicalc.fxi,K,Kvals);
    u0 = replace(ulaw.u0,K,Kvals);
    ucomp = replace(ulaw.ucomp,K,Kvals);


    if isempty(Qvals)           
        V = inicalc.V; 
        dVdx = inicalc.dVdx;
    else   
        V = mV.'*Qvals*mV; 
        dVdx = clean(jacobian(V,xi), numsets.clean_thresh);        
    end

    dVdt = dVdx*fxi;

    prob = [];
    %expand conditions
    eq1 = - dVdt  - s1 * (alpha-V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2  + (V - alpha))- numsets.epsi;
    F1 = [sos(s1),sos(eq1)];
    F2 = [sos(s2),sos(eq2)];
    [sol] = solvesos(F1,[],numsets.sdpsetting1,c1);
    prob = [prob,sol.problem];
    [sol] = solvesos(F2,[],numsets.sdpsetting1,c2);
    prob = [prob,sol.problem];
    if any(prob>0) ||any(prob<-1)
        fprintf(1,'The following errors were found: %d %d \n',[prob]);
        error('Error in multiplier search')
    elseif any(prob==-1)
        if length(Smatlist)>=5
            multis = multisold;
            multis.flag =-1;
            fprintf(1,'-------Ran into numerical problems, final results file created but terminating early..-------\n\n')
            return
        else
            fprintf(1,'The following errors were found: %d %d \n',[prob]);
            error('Error in multiplier search')
        end
    else
        multis.c1 = value(c1)'; 
        multis.c2 = value(c2)';   
    end
    
    
    for i = 1:sys.udim
       [s3,c3,~] = polynomial([xi],degs.s3_dU,degs.s3_dL);
       [s4,c4,~] = polynomial([xi],degs.s4_dU,degs.s4_dL);
        s5 = 1;
        s6 = 1;
        eq3a = -s3*(u0(i) + sum(tvar(i,:)) - sys.umax) - s5*(alpha-V) - numsets.epsi;
        eq3b = -s4*(sys.umin - (u0(i)-(sum(tvar(i,:))))) - s6*(alpha-V) - numsets.epsi;
        F3 = [sos(s3),sos(eq3a)]; 
        [sol] = solvesos(F3,[],numsets.sdpsetting1,[c3(:)]);         
        prob = [prob,sol.problem];
        F4 = [sos(s4),sos(eq3b)];
        [sol] = solvesos(F4,[],numsets.sdpsetting1,[c4(:)]);
        prob = [prob,sol.problem];
        
        if any(prob>0) ||any(prob<-1)
            fprintf(1,'The following errors were found: %d %d \n',[prob]);
            error('Error in multiplier search')
        elseif any(prob==-1) 
            if length(Smatlist)>=5
                multis = multisold;
                multis.flag =-1;
                fprintf(1,'-------Ran into numerical problems, final results file created but terminating early..-------\n\n')
                return
            else
                fprintf(1,'The following errors were found: %d %d \n',[prob]);
                error('Error in multiplier search')
            end
        else
            multis.c3{i} = value(c3)'; 
            multis.c4{i} = value(c4)';   
        end  
        
        for j = 1:(sys.p)
            [s7,c7,~] = polynomial([xi],degs.s7_dU,degs.s7_dL);
            [s8,c8,~] = polynomial([xi],degs.s8_dU,degs.s8_dL);
            s9 = 1;
            s10 = 1;
            eq4a = s7*(tvar(i,j)+ucomp(i,j)) - s9*(alpha-V) - numsets.epsi;
            eq4b = s8*(tvar(i,j)-ucomp(i,j)) - s10*(alpha-V) - numsets.epsi;
            F5 = [sos(s7),sos(eq4a)]; 
            [sol] = solvesos(F5,[],numsets.sdpsetting1,[c7(:)]);
            prob = [prob,sol.problem];
            F6 = [sos(s8),sos(eq4b)];
            [sol] = solvesos(F6,[],numsets.sdpsetting1,[c8(:)]); 
            prob = [prob,sol.problem];
            if any(prob>0) ||any(prob<-1)
                fprintf(1,'The following errors were found: %d %d \n',[prob]);
                error('Error in multiplier search')
            elseif any(prob==-1)
                if length(Smatlist)>=5
                    multis = multisold;
                    multis.flag =-1;
                    fprintf(1,'-------Ran into numerical problems, final results file created but terminating early..-------\n\n')
                    return
                else
                    fprintf(1,'The following errors were found: %d %d \n',[prob]);
                    error('Error in multiplier search')
                end
            else
                multis.c7{i}{j} = value(c7)'; 
                multis.c8{i}{j} = value(c8)';  
            end
        
        end
    end
            
    multis.flag = 0;                 
    disp('-------Multiplier search successfully finished-------')


end


function [Smatres,Qvalres] = optimize_V_wIC(xi,K,sys,numsets,degs,inicalc,multis,Smatold,gterm,Smatlist,Qvalsold,Kvals,tvar,ulaw)

    gamma = 1;
    alpha = 1;
    
    fxi = replace(inicalc.fxi,K,Kvals);
    u0 = replace(ulaw.u0,K,Kvals);
    ucomp = replace(ulaw.ucomp,K,Kvals);
        
    Smatdiff = sdpvar(length(xi),length(xi));
    Smat = Smatold - Smatdiff;
  
    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);       
    
    [~,~,m1] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([xi],degs.s2_dU,degs.s2_dL);
    [~,~,m3] = polynomial([xi],degs.s3_dU,degs.s3_dL);
    [~,~,m4] = polynomial([xi],degs.s4_dU,degs.s4_dL);
    [~,~,m7] = polynomial([xi],degs.s7_dU,degs.s7_dL);
    [~,~,m8] = polynomial([xi],degs.s8_dU,degs.s8_dL);


   
    Qvals = sdpvar(size(mV,1),size(mV,1)); 
    V = mV.'*Qvals*mV;      
    dVdx = clean(jacobian(V,xi), numsets.clean_thresh);
    dVdt = dVdx*fxi;
    
    F = [];
    s1fix = multis.c1*m1;
    s2fix = multis.c2*m2;
    
    eq1 = - dVdt - s1fix * (alpha - V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2fix + (V - alpha))- numsets.epsi; 
    F = [F,sos(eq1),sos(eq2)];    
    
    for i = 1:sys.udim
        s3fix = multis.c3{i}*m3;
        s4fix = multis.c4{i}*m4;
        s5fix = 1;
        s6fix = 1;
        eq3a = -s3fix*(u0(i) + sum(tvar(i,:)) - sys.umax) - s5fix*(alpha-V) - numsets.epsi;
        eq3b = -s4fix*(sys.umin - (u0(i)-(sum(tvar(i,:))))) - s6fix*(alpha-V) - numsets.epsi;
        F = [F,sos(eq3a),sos(eq3b)]; 
        
        for j = 1:(sys.p)
            s7fix = multis.c7{i}{j}*m7;
            s8fix = multis.c8{i}{j}*m8;
            s9fix = 1;
            s10fix = 1;
            eq4a = s7fix*(tvar(i,j)+ucomp(i,j)) - s9fix*(alpha-V) - numsets.epsi;
            eq4b = s8fix*(tvar(i,j)-ucomp(i,j)) - s10fix*(alpha-V) - numsets.epsi;
            F = [F,sos(eq4a),sos(eq4b)]; 
        end
    end
    
    F = [F,Qvals>=0,Smat>=0,Smatdiff>=0];
 
    [sol] = solvesos(F, -geomean(Smatdiff), numsets.sdpsetting2,[Qvals(:);Smatdiff(:)]);
    
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





function [Kvalres,tvarres] = optimize_K_wIC(xi,K,sys,numsets,degs,inicalc,multis,Smat,gterm,Qvals,ulaw)
    
    gamma = 1;
    alpha = 1;
    
    tvar = sdpvar(sys.udim,sys.p);
 
    fxi = inicalc.fxi;
    u0 = ulaw.u0;
    ucomp = ulaw.ucomp;
 
    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);       
    
    [~,~,m1] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([xi],degs.s2_dU,degs.s2_dL);
    [~,~,m3] = polynomial([xi],degs.s3_dU,degs.s3_dL);
    [~,~,m4] = polynomial([xi],degs.s4_dU,degs.s4_dL);
    [~,~,m7] = polynomial([xi],degs.s7_dU,degs.s7_dL);
    [~,~,m8] = polynomial([xi],degs.s8_dU,degs.s8_dL);
  
    if isempty(Qvals)           
        error('Qvals cannot be empty at this point!')
    else   
        V = mV.'*Qvals*mV; 
        dVdx = clean(jacobian(V,xi), numsets.clean_thresh);        
    end

    dVdt = dVdx*fxi;
    
    F = [];
    s1fix = multis.c1*m1;
    s2fix = multis.c2*m2;
    
    eq1 = - dVdt - s1fix * (alpha - V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2fix + (V - alpha))- numsets.epsi; 
    F = [F,sos(eq1),sos(eq2)];    
    
    for i = 1:sys.udim
        s3fix = multis.c3{i}*m3;
        s4fix = multis.c4{i}*m4;
        s5fix = 1;
        s6fix = 1;
        eq3a = -s3fix*(u0(i) + sum(tvar(i,:)) - sys.umax) - s5fix*(alpha-V) - numsets.epsi;
        eq3b = -s4fix*(sys.umin - (u0(i)-(sum(tvar(i,:))))) - s6fix*(alpha-V) - numsets.epsi;
        F = [F,sos(eq3a),sos(eq3b)]; 
        
        for j = 1:(sys.p)
            s7fix = multis.c7{i}{j}*m7;
            s8fix = multis.c8{i}{j}*m8;
            s9fix = 1;
            s10fix = 1;
            eq4a = s7fix*(tvar(i,j)+ucomp(i,j)) - s9fix*(alpha-V) - numsets.epsi;
            eq4b = s8fix*(tvar(i,j)-ucomp(i,j)) - s10fix*(alpha-V) - numsets.epsi;
            F = [F,sos(eq4a),sos(eq4b),tvar(i,j)>=0]; 
        end
    end
    
 
    [sol] = solvesos(F, [], numsets.sdpsettingK,[K(:);tvar(:)]);
    
    if sol.problem ~= 0 
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in K optimization')
    end
    
    Kvalres = value(K);
    tvarres = value(tvar);
     
    fprintf(1,'-------K optimization successfully finished-------\n\n')
     

end



function [Kvalres] = optimize_K(xi,K,numsets,degs,Qvals,inicalc,gterm,multis)
  
    alpha = 1;

    [~,~,m1] = polynomial([xi],degs.s1_dU,degs.s1_dL);

    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);
    
    if isempty(Qvals)           
        V = inicalc.V; 
        dVdx = inicalc.dVdx;
    else   
        V = mV.'*Qvals*mV; 
        dVdx = clean(jacobian(V,xi), numsets.clean_thresh);        
    end
       
    Vdot = dVdx*inicalc.fxi;
    s1fix = multis.c1*m1;

    eq1 = - Vdot - s1fix * (alpha - V) - gterm;

    F = [sos(eq1)];
       
    sol = solvesos(F, [], numsets.sdpsettingK,[K(:)]);
    
    if sol.problem ~= 0 
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in K optimization')
    end
    
    Kvalres = value(K);
     
    fprintf(1,'-------K optimization successfully finished-------\n\n')
        

end

function multis = find_multiplier(xi,K,numsets,degs,Qvals,inicalc,gterm,Smat,Smatlist,multisold,Kvals)

    alpha = 1;
    gamma = 1;
                 
    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);
    
    [s1,c1,~] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [s2,c2,~] = polynomial([xi],degs.s2_dU,degs.s2_dL);
    
    fxi = replace(inicalc.fxi,K,Kvals);
    
 
    if isempty(Qvals)           
        V = inicalc.V; 
        dVdx = inicalc.dVdx;
    else   
        V = mV.'*Qvals*mV; 
        dVdx = clean(jacobian(V,xi), numsets.clean_thresh);        
    end
    
    dVdt = dVdx*fxi;
    
    eq1 = - dVdt - s1 * (alpha-V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2 + (V - alpha))- numsets.epsi;
           
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
        disp([res1,res2]);
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


function [Smatres,Qvalres] = optimize_V(xi,K,numsets,degs,inicalc,multis,Smatold,gterm,Smatlist,Qvalsold,Kvals)

    alpha = 1;
    gamma = 1;
    
    fxi = replace(inicalc.fxi,K,Kvals);
    
    Smatdiff = sdpvar(length(xi),length(xi));
    Smat = Smatold - Smatdiff;
    
    [~,~,m1] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([xi],degs.s2_dU,degs.s2_dL);
   
    mV = monolist(xi,degs.V_dU/2,degs.V_dL/2);
       
    Qvals = sdpvar(size(mV,1),size(mV,1)); 
    V = mV.'*Qvals*mV;      
    dVdx = clean(jacobian(V,xi), numsets.clean_thresh);
    dVdt = dVdx*fxi;
               
    s1fix = multis.c1*m1;
    s2fix = multis.c2*m2;
    
    eq1 = - dVdt - s1fix * (alpha - V) - gterm ;
    eq2 = -((gamma - xi.'*Smat*xi) * s2fix + (V - alpha))- numsets.epsi; 

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






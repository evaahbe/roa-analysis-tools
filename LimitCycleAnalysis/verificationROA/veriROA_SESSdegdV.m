function veriROA_SESSdegdV(sys,fns,filenames,pinit,numsets)

    
    rho = sdpvar(sys.xdim-1,1); % create n-1 dim state on hyperplane
    
    tau_array = linspace(0,sys.Tperiod,sys.Nhp)'; % create array of tau samples
    
    degs = numsets.degs;

    % load or compute solution of periodic Lyapunov/Riccati equation
    try
        load(filenames.periodinitname)
    catch
        if sys.udim ==0
            [Pvals, Pdotvals] = fns.periodinitgain(tau_array,sys,fns,pinit);
            save(filenames.periodinitname, 'Pvals', 'Pdotvals')
        else
            [Pvals, Pdotvals, Kvals] = fns.periodinitgain(tau_array,sys,fns,pinit);
            save(filenames.periodinitname, 'Pvals', 'Pdotvals', 'Kvals')
        end
    end
    
    % compute the transverse dynamics on each hyperplane 
    for i = 1:sys.Nhp 
        tau = tau_array(i);
        
        xorbi = ppval(sys.xorb, tau);
        [Zi, Zdoti] = projectOperator(sys, tau);
        xi = xorbi + Zi' * rho; 
        
        if sys.udim == 0 
            % scale solution of Lyapunov equation such that initially feasible
            Pmati = Pvals{i}/pinit.initVscale;
            Pmatdoti = Pdotvals{i}/pinit.initVscale;
            
            forbi = fns.system_dyn(xorbi,sys);
            [fxi] = fns.system_dyn_poly(xi,sys);     
        else 
            % scale solution of Riccati equation such that initially feasible
            Pmati = Pvals{i}/pinit.initVscale;
            Pmatdoti = Pdotvals{i}/pinit.initVscale;
            Ki = Kvals{i};
            
            uorbi = ppval(sys.uorb, tau);
            forbi = fns.system_dyn(xorbi,uorbi,sys);
            [fxi] = fns.system_dyn_poly(xorbi,uorbi,sys,xi,Ki,Zi); 
        end
        vi = ppval(sys.v,tau);
        
        taudotdeni = vi'*forbi + vi'*Zdoti'*rho;
        taudotdeni = clean(taudotdeni,numsets.clean_thresh);
        taudotnumi = vi'*fxi;
        taudotnumi = clean(taudotnumi,numsets.clean_thresh);
        
        rhodoti = clean(Zdoti*Zi'*rho*taudotnumi,numsets.clean_thresh) + clean(Zi*fxi*taudotdeni,numsets.clean_thresh) - clean(Zi*forbi*taudotnumi,numsets.clean_thresh); 
        rhodoti = clean(rhodoti,numsets.clean_thresh);

        % save in struct for direct access in iteration 
        onHP{i}.rhodot      = rhodoti;
        onHP{i}.tau         = tau;
        onHP{i}.taudotden   = taudotdeni;
        onHP{i}.taudotnum   = taudotnumi;
        onHP{i}.Qini        = Pmati;
        onHP{i}.Qinidot     = Pmatdoti;

    end
    
    try             % load results of already performed iteration steps
        load(filenames.intermedresults)
    catch           % create intermediate results file for possible reinitialization
        %initializing the variables
        step    = 1;                % current step, for stopping and reinitializing iteration 
        Qvals   = [];               % Gram matrix entries of V
        multis  = [];               % coefficients of SOS and other multiplier
        gamma   = 1;                % V sublevel set, in SESSdegdV this is fixed
        alpha   = numsets.alpha;    % B sublevel set, in SESSdegdV this is optimized over
        alpha_list = [-alpha,alpha]; % initilize alpha list 
        save(filenames.intermedresults,'alpha_list','step','Qvals','alpha','gamma','multis')
    end
    
    
    lterm = numsets.lfac * rho'*rho;
     
    jj=length(alpha_list(:))-1;
    while jj < numsets.iteration_max
        
        %check if convergence criteria satisfied
        if (alpha_list(end) - alpha_list(end-1)) < numsets.convCrit*alpha_list(end-1)         
            disp('Convergence reached!')
            break;
        end
        
        fprintf(1,'\n \n This is iteration nr.: %d \n',jj)
        
        if step == 1                % feasibility test for multiplier
            fprintf(1,'--> multiplier test \n')
            multis = find_multiplier(rho,onHP,sys,numsets,degs,lterm,Qvals,alpha);
            step = 2;
            save(filenames.intermedresults,'-append','step','multis')
        elseif step == 2            % maximize ROA volume by searching over V
            fprintf(1,'--> V optimization \n')
            [alpha,Qvals] = optimize_V(rho,onHP,sys,numsets,degs,lterm,multis,alpha_list);
            alpha_list = [alpha_list,alpha];
            step = 1;
            jj = jj+1;
            save(filenames.intermedresults,'-append','step','alpha','Qvals','alpha_list')
        end  

    end
    save(filenames.finalresults, 'multis','Qvals','step','alpha_list','gamma','alpha','sys','numsets')  
end


function multis = find_multiplier(rho,onHP,sys,numsets,degs,lterm,Qvals,alpha)

    
    Smatfix = numsets.ellipsemat;
    
    mV = monolist(rho,degs.V_dU/2,degs.V_dL/2);
    
    for i = 1:sys.Nhp-1
        
        % initilize SOS multiplier
        [s1,c1, ~] = polynomial([rho],degs.s1_dU,degs.s1_dL);
        [s2,c2, ~] = polynomial([rho],degs.s2_dU,degs.s2_dL);
        [s3,c3, ~] = polynomial([rho],degs.s3_dU,degs.s3_dL);

        prob = [];      
        % rebuild fixed V
        if isempty(Qvals)   % in the very first iteration use Lyapunov/Riccati solution

            V = rho' * onHP{i}.Qini * rho;
            delVdelrho = clean(jacobian(V,rho), numsets.clean_thresh); 

            delVdeltau = rho' * onHP{i}.Qinidot * rho;
        else                % then, use result from step 2
            if i~=sys.Nhp-1
                Qmat = reshape(Qvals(i,:),size(mV,1),size(mV,1));               
                V = mV'*Qmat*mV;

                Qmatp = reshape(Qvals(i+1,:),size(mV,1),size(mV,1));
                Vp = mV'*Qmatp*mV;

                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;

                delVdeltau = (Vp-V)/(taup-tau);
                if i==1
                    VN = V;
                end     
            else            % for correct (periodic) end point derivative 
                Qmat = reshape(Qvals(i,:),size(mV,1),size(mV,1));
                V = mV'*Qmat*mV;
                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;

                delVdeltau = (VN-V)/(taup-tau);        
            end   
            delVdelrho = clean(jacobian(V,rho), numsets.clean_thresh);
        end
        
        % total Lyapunov derivative
        dVdt = delVdelrho*onHP{i}.rhodot + delVdeltau*onHP{i}.taudotnum + numsets.neg_def_fac*rho'*rho*onHP{i}.taudotden; 
        
        %%% SOS PROGRAM
        % ROA conditions on a single hyperplane 
        eq1 = - dVdt - s1 * (1-V)-lterm ;    
        eq2 = onHP{i}.taudotden - s2 * (1-V)- numsets.epsi;
        eq3 = -((alpha - rho'*Smatfix*rho) * s3 + (V - 1));

        F1 = [sos(eq1),sos(s1)];
        F2 = [sos(eq2),sos(s2)];
        F3 = [sos(eq3),sos(s3)];

        % solve for each constraint separately to obtain multiplier   
        [sol1] = solvesos(F1, [], numsets.sdpsetting1, c1); 
        [sol2] = solvesos(F2, [], numsets.sdpsetting1, c2);
        [sol3] = solvesos(F3, [], numsets.sdpsetting1, c3);

        prob = [prob,sol1.problem];
        prob = [prob,sol2.problem];
        prob = [prob,sol3.problem];
        

        % in case feasibility test is infeasible 
        if any(prob~=0)
            fprintf(1,'In round %d the following errors were found: %d %d %d\n',[i,prob]);
            error('Error in multiplier search')
        end
        
        % collect results to return
        multis{i}.c1 = value(c1)'; 
        multis{i}.c2 = value(c2)';
        multis{i}.c3 = value(c3)';
    end
                      
    fprintf(1,'-------Multiplier search successfully finished-------\n\n')
end


function [alphares,Qvals] = optimize_V(rho,onHP,sys,numsets,degs,lterm,multis,alpha_list)


    Smatfix = numsets.ellipsemat;
    
    [~,~,m1] = polynomial([rho],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([rho],degs.s2_dU,degs.s2_dL);
    [~,~,m3] = polynomial([rho],degs.s3_dU,degs.s3_dL);
     
    F  = []; 
    
    mV = monolist(rho,degs.V_dU/2,degs.V_dL/2);
    Qcols = size(mV,1)^2;
    
    alpha = sdpvar(1,1);
    alphaold = alpha_list(end);
    alphatot = alphaold + alpha;
        
    for i=1:sys.Nhp-1
      
        % rebuild fixed multiplier
        s1fix = clean(multis{i}.c1*m1,numsets.clean_thresh);
        s2fix = clean(multis{i}.c2*m2,numsets.clean_thresh);
        s3fix = clean(multis{i}.c3*m3,numsets.clean_thresh);

        % initialize V
        if i < sys.Nhp-1
            if i ==1        % for periodicity
                Qmat = sdpvar(size(mV,1),size(mV,1));
                Qmat1 = Qmat;
                F = [F,Qmat>=0];
                Qfull(i,:) = reshape(Qmat,1,Qcols);
            else
                Qmat = Qmatp;
            end
            Qmatp = sdpvar(size(mV,1),size(mV,1)); 

            V = mV'*Qmat*mV;      
            Vp = mV'*Qmatp*mV;

            tau = onHP{i}.tau;
            taup = onHP{i+1}.tau;

            delVdeltau = (Vp-V)/(taup-tau);
            delVdelrho = clean(jacobian(V,rho), numsets.clean_thresh);
            
            Qfull(i+1,:) = reshape(Qmatp,1,Qcols);
        else                % make sure periodicity is maintained     
            Qmat = Qmatp; 
            V = mV'*Qmat*mV;
            Vp = mV'*Qmat1*mV; 
            tau = onHP{i}.tau;
            taup = onHP{i+1}.tau;

            delVdeltau = (Vp-V)/(taup-tau);
            delVdelrho = clean(jacobian(V,rho), numsets.clean_thresh);
        end     
        
        dVdt = delVdelrho*onHP{i}.rhodot + delVdeltau*onHP{i}.taudotnum + numsets.neg_def_fac*rho'*rho*onHP{i}.taudotden; 
        
        %SOS PROGRAM
        % collect ROA conditions for each hyperplane
        eq1 = - dVdt - s1fix * (1 - V) -lterm ;
        eq2 = onHP{i}.taudotden - s2fix * (1 - V)- numsets.epsi;
        eq3 = -((alphatot - rho'*Smatfix*rho) * s3fix + (V - 1)); 

        F = [F, sos(eq1), sos(eq2), sos(eq3),Qmatp>=0];
        
    end
    
    F = [F,alpha>=0];

    % solve SOS Program
    [sol] = solvesos(F, -alpha, numsets.sdpsetting2 ,[Qfull(:);alpha]);
    
    if sol.problem ~= 0
        disp(sol.problem)
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in V optimization')
    end
    
    % rebuild results to return
    Qvals = value(Qfull);
    alphares = value(alphatot);

    fprintf('This is the alpha result: %f \n',alphares);

    fprintf(1,'-------V optimization successfully finished-------\n\n')
end







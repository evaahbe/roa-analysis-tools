function veriROA_EEdegdV(sys,fns,filenames,pinit,numsets)

    
    rho = sdpvar(sys.xdim-1,1); % create n-1 dim state on hyperplane
    
    tau_array = linspace(0,sys.Tperiod,sys.Nhp)'; % create array of tau samples
    
    degs = numsets.degs;

    % load or compute solution of periodic Lyapunov/Riccati equation
    try
        load(filenames.periodinitname)
    catch
        if sys.udim ==0
            %[Pvals, Pdotvals] = devalPeriodicLyapunov_testing(tau_array,sys,fns,pinit);
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
        % initialize the surrogate set sizes on the hyperplanes
        [Smat] = createInitialSmat(onHP,sys); 
        if numsets.tvball == 1          %take surrogate varying with hyperplanes
            Smatlist(:,:,1) = 10*repmat(reshape(Smat,1,length(rho)^2),sys.Nhp-1,1);
            Smatlist(:,:,2) = 1/10*Smatlist(:,:,1);
            geomeanS_list(:,1) = 10*repmat(det(Smat)^(1/length(Smat)),sys.Nhp-1,1);
            geomeanS_list(:,2) = 1/10* geomeanS_list(:,1);
        elseif numsets.tvball == 0      %take the same surrogate for all hyperplanes
            Smatlist(1,:) = 10*reshape(Smat,1,length(rho)^2);
            Smatlist(2,:) = 1/10*Smatlist(1,:);
            geomeanS_list(1) = 10*det(Smat)^(1/length(Smat));
            geomeanS_list(2) = 1/10* geomeanS_list(1);
        end
        %initializing the variables
        step    = 1;  % current step, for stopping and reinitializing iteration 
        Qvals   = []; % Gram matrix entries of V
        multis  = []; % coefficients of SOS and other multiplier
        gamma   = 1;  % V sublevel set, in EEdegdV this is fixed
        alpha   = 1;  % B sublevel set, in EEdegdV this is fixed
        save(filenames.intermedresults,'Smatlist','geomeanS_list','step','Qvals','gamma','alpha','multis')
    end
       
    lterm = numsets.lfac * rho'*rho;
    
    if numsets.tvball == 1
        jj=length(geomeanS_list(1,:))-1;
    elseif numsets.tvball == 0
        jj=length(geomeanS_list(:))-1;
    end
    while jj < numsets.iteration_max
        
        %check if convergence criteria satisfied
        if numsets.tvball ==1
            if (max(geomeanS_list(:,end-1) - geomeanS_list(:,end))) < numsets.convCrit*max(geomeanS_list(:,end))  
                disp('Convergence reached!')
                break;
            end
        elseif numsets.tvball ==0
            if (geomeanS_list(end-1) - geomeanS_list(end)) < numsets.convCrit*geomeanS_list(end) 
                disp('Convergence reached!')
                break;
            end  
        end
            
        fprintf(1,'\n \n This is iteration nr.: %d \n',jj)
        
        if step == 1                % feasibility test for multiplier
            fprintf(1,'--> multiplier test \n')
            multis = find_multiplier(rho,onHP,sys,numsets,degs,lterm,Qvals,Smatlist);
            step = 2;
            save(filenames.intermedresults,'-append','step','multis')
        elseif step == 2            % maximize ROA volume by searching over V
            fprintf(1,'--> V optimization \n')
            [Smat,geomean_Smat,Qvals] = optimize_V(rho,onHP,sys,numsets,degs,lterm,multis,Smatlist);
            if numsets.tvball ==1
                Smatlist(:,:,end+1) = Smat;
                geomeanS_list(:,size(geomeanS_list,2)+1) = geomean_Smat(:);
            elseif numsets.tvball ==0
                Smatlist = [Smatlist;Smat];
                geomeanS_list = [geomeanS_list,geomean_Smat];
            end
            step = 1;
            jj = jj+1;
            save(filenames.intermedresults,'-append','step','Smatlist','Qvals','geomeanS_list')
        end  

    end
    save(filenames.finalresults, 'multis','Smatlist','Qvals','step','geomeanS_list','gamma','alpha','sys','numsets')  
end


function multis = find_multiplier(rho,onHP,sys,numsets,degs,lterm,Qvals,Smatlist)

    alpha = 1;
    gamma = 1;
    
    mV = monolist(rho,degs.V_dU/2,degs.V_dL/2);
    
    for i = 1:sys.Nhp-1
        
        % rebuild fixed surrogate set 
        if numsets.tvball == 1
            Smat = reshape(Smatlist(i,:,end),length(rho),length(rho));
        elseif numsets.tvball == 0
            Smat = reshape(Smatlist(end,:),length(rho),length(rho));
        end
        
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
        dVdt = delVdelrho*onHP{i}.rhodot + delVdeltau*onHP{i}.taudotnum + numsets.neg_def_fac*(rho'*rho)*onHP{i}.taudotden; 
        
        %%% SOS PROGRAM
        % ROA conditions on a single hyperplane 
        eq1 = - dVdt - s1 * (gamma-V)-lterm ;    
        eq2 = onHP{i}.taudotden - s2 * (gamma-V)- numsets.epsi;
        eq3 = -((alpha - rho'*Smat*rho) * s3 + (V - gamma));

        F1 = [sos(eq1),sos(s1)];
        F2 = [sos(eq2),sos(s2)];
        F3 = [sos(eq3),sos(s3)];

        % solve for each constraint separately to obtain multiplier   
        [sol1] = solvesos(F1, [], numsets.sdpsetting1, c1); 
        [sol2] = solvesos(F2, [], numsets.sdpsetting1, c2);
        [sol3] = solvesos(F3, [], numsets.sdpsetting1, c3);

        res = sol1.problem;
        prob = [prob,res];
        res = sol2.problem;
        prob = [prob,res];
        res = sol3.problem;
        prob = [prob,res];
        

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


function [Smatvals,geomean_Smat,Qvals] = optimize_V(rho,onHP,sys,numsets,degs,lterm,multis,Smatlist)


    gamma = 1;
    alpha = 1;
    
    [~,~,m1] = polynomial([rho],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([rho],degs.s2_dU,degs.s2_dL);
    [~,~,m3] = polynomial([rho],degs.s3_dU,degs.s3_dL);
     
    F           = []; 
    ObjFun      = [];
    
    mV = monolist(rho,degs.V_dU/2,degs.V_dL/2);
    Qcols = size(mV,1)^2;
        
    for i=1:sys.Nhp-1
   
        % initialize surrogate set for optimization (to maximize the
        % DIFFERENCE from the previous size)
        if numsets.tvball == 1 
            Smat = sdpvar(length(rho),length(rho)); 
            Smatold = reshape(Smatlist(i,:,end),length(rho),length(rho));
            Smattot = Smatold - Smat;
            Smatfull(i,:) = reshape(Smat,1,length(rho)^2);
            Smatoldlist(i,:) = Smatlist(i,:,end);
            
            F = [F,Smat>=0,Smattot>=0];
            ObjFun = [ObjFun,-geomean(Smat)];
            
        elseif numsets.tvball == 0 
            if i == 1
                Smat = sdpvar(length(rho),length(rho));
                Smatold = reshape(Smatlist(end,:),length(rho),length(rho));
                Smattot = Smatold - Smat;
                Smatoldlist = Smatlist(end,:);
                Smatfull = reshape(Smat,1,length(rho)^2);
                
                F = [F,Smat>=0,Smattot>=0];
                ObjFun = [ObjFun,-geomean(Smat)];
                
            end
        end       
        
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
        
        dVdt = delVdelrho*onHP{i}.rhodot + delVdeltau*onHP{i}.taudotnum + numsets.neg_def_fac*(rho'*rho)*onHP{i}.taudotden; 
        
        %SOS PROGRAM
        % collect ROA conditions for each hyperplane
        eq1 = - dVdt - s1fix * (gamma - V) -lterm ;
        eq2 = onHP{i}.taudotden - s2fix * (gamma - V)- numsets.epsi;
        eq3 = -((alpha - rho'*Smattot*rho) * s3fix + (V - gamma)); 

        F = [F, sos(eq1), sos(eq2), sos(eq3),Qmatp>=0];
        
    end
  
    % solve SOS Program
    [sol] = solvesos(F, ObjFun, numsets.sdpsetting2,[Qfull(:);Smatfull(:)]);
    
    if sol.problem ~= 0
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in V optimization')
    end
    
    % rebuild results to return
    Qvals = value(Qfull);
    Smatvals = Smatoldlist - value(Smatfull);
    
    if numsets.tvball ==1
        for i = 1:sys.Nhp-1
            Smatres = reshape(Smatvals(i,:),length(rho),length(rho));
            geomean_Smat(i) = det(Smatres)^(1/length(Smatres)); %formula for geomean
        end
    elseif numsets.tvball ==0
        Smatres = reshape(Smatvals,length(rho),length(rho));
        geomean_Smat = det(Smatres)^(1/length(Smatres)); %formula for geomean
        fprintf(1,'This is the geomean result: %f \n',geomean_Smat);
    end
    
    fprintf(1,'-------V optimization successfully finished-------\n\n')
end




function [Smat] = createInitialSmat(onHP,sys)
    
    % creates a feasible initial surrogate size to fit inside initial V.
    for i = 1:length(onHP)
        eigenvals(i) = max(eig(onHP{i}.Qini));
    end
    
    Smat = eye(sys.xdim-1)*max(eigenvals)*1.2; 
end





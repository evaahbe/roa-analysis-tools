function veriROC_maxROCballInv(sys,fns,filenames,pinit,numsets,~)

    rho = sdpvar(sys.xdim-1,1); % create n-1 dim state on hyperplane
    y = sdpvar(sys.xdim-1,1); % create n-1 dim indeterminant vector for matrix SOS check
    
    tau_array = linspace(0,sys.Tperiod,sys.Nhp)'; % create array of tau samples
    
    degs = numsets.degs;
    
    % load or compute solution of periodic Lyapunov/Riccati equation
    try
        load(filenames.periodinitname)
    catch
        if sys.udim ==0
            [Pvals, Pdotvals] = fns.periodinitgain(tau_array,sys,fns,pinit);
            save(filenames.periodinitname, 'Pvals', 'Pdotvals')
            Kvals = [];
        else
            [Pvals, Pdotvals, Kvals] = fns.periodinitgain(tau_array,sys,fns,pinit);
            save(filenames.periodinitname, 'Pvals', 'Pdotvals', 'Kvals')
        end
    end
        
 
    for i = 1:sys.Nhp
        tau = tau_array(i);
             
        xorbi = ppval(sys.xorb, tau);
        [Zi, Zdoti] = projectOperator(sys, tau);
        xloci = xorbi + Zi' * rho; 
        vi = ppval(sys.v,tau);
        
        if sys.udim == 0 
            % scale solution of Lyapunov equation such that initially feasible
            Pmati = Pvals{i}/pinit.initVscale;
            Pmatdoti = Pdotvals{i}/pinit.initVscale;
            
            forbi = fns.system_dyn(xorbi,sys);  
            
            onHP{i}.K = 0;
            Kvals = [];
        else 
            % scale solution of Riccati equation such that initially feasible
            Pmati = Pvals{i}/pinit.initVscale;
            Pmatdoti = Pdotvals{i}/pinit.initVscale;
            
            uorbi = ppval(sys.uorb, tau);
            forbi = fns.system_dyn(xorbi,uorbi,sys);

            onHP{i}.uorb = uorbi;
            onHP{i}.K = Kvals{i};
        end
       

        taudotdeni = vi'*forbi + vi'*Zdoti'*rho;
        taudotdeni = clean(taudotdeni,numsets.clean_thresh);

        % save in struct for direct access in iteration 
        onHP{i}.tau         = tau;
        onHP{i}.taudotden   = taudotdeni;
        onHP{i}.Qini        = Pmati;
        onHP{i}.Qinidot     = Pmatdoti;
        onHP{i}.xloc        = xloci;
        onHP{i}.forb        = forbi;
        onHP{i}.xorb        = xorbi;
        onHP{i}.Z           = Zi;
        onHP{i}.Zdot        = Zdoti;
        onHP{i}.v           = vi;

    end
    
    Smatfix = numsets.ellipsemat; 
      

    %initializing the variables
    step        = 1;                % current step, for stopping and reinitializing iteration 
    Qvals       = [];               % Gram matrix entries of V
    Bvals       = [];               % Gram matrix entries of B
    multis      = [];               % coefficients of SOS and other multiplier
    multi_mats  = [];               % coefficients of the matrix SOS multiplier
    gamma       = 1;                % V sublevel set, fixed to a constant here 
    alpha       = numsets.alpha;    % B sublevel set, in maxROC this is optimized for

    if numsets.tvalpha == 1         %take sublevel set constant varying with hyperplanes
        alpha = repmat(alpha,1,sys.Nhp-1);
        alpha_list(1,:) = zeros(size(alpha));
        alpha_list(2,:) = alpha;
    elseif numsets.tvalpha == 0     %take the same sublevel set constant for all hyperplanes
        alpha_list(1,1) = 0;
        alpha_list(2,1) = alpha;
    end

    try             % load results of already performed iteration steps
        load(filenames.rocintermedresults)
    catch           % create intermediate results file for possible reinitialization
        save(filenames.rocintermedresults,'alpha_list','step','Qvals','gamma','multis','multi_mats','alpha')
    end
    
    jj = length(alpha_list(:,1))-1;
    while jj<=numsets.iteration_max
        
        if (max(alpha_list(end,:) - alpha_list(end-1,:))) < numsets.convCrit*max(alpha_list(end-1,:))
            disp('Convergence reached!')
            break;
        end
        
        fprintf(1,'\n \n This is iteration nr.: %d \n',jj);
        if step == 1   
            fprintf(1,'- step: Multiplier search\n\n')
            [multis, multi_mats] = find_multiplier(y,rho,fns,onHP,sys,numsets,degs,Qvals,alpha,Bvals);
            step = 2;
            save(filenames.rocintermedresults,'-append','step','multis','multi_mats')
        elseif step == 2 
            fprintf(1,'- step: Metric optimization\n\n')
            [alpha,Qvals,Bvals] = optimize_Metric(y,rho,fns,onHP,sys, numsets,degs,alpha_list,multis,multi_mats);
            alpha_list = [alpha_list;alpha];
            step = 1;
            jj = jj+1;
            save(filenames.rocintermedresults,'-append','alpha_list','Qvals','step','alpha','Bvals')
        end
            
    end
    Qvalsmet = Qvals;
    Qvals = Bvals;
    save(filenames.rocfinalresults, 'multis','multi_mats','alpha_list','Qvals','alpha','Kvals','gamma','sys','numsets','Qvalsmet')   
end



function [multis,multi_mats] = find_multiplier(y,rho,fns,onHP,sys,numsets,degs,Qvals,alpha,Bvals)


    Smatfix = numsets.ellipsemat;
    [~,~,mM] = polynomial(rho,degs.M_dU,degs.M_dL);
    
    for i = 1:sys.Nhp-1
        prob = [];
        
        taudotden = onHP{i}.taudotden;
        
        [s2,c2, ~] = polynomial(rho,degs.s2_dU,degs.s2_dL);
        [s3,c3, ~] = polynomial(rho,degs.s3_dU,degs.s3_dL);
        
        if numsets.tvalpha == 1          
            alphai = alpha(i);
        elseif numsets.tvalpha == 0      
            alphai = alpha;
        end
      
        if isempty(Qvals)
            
            M = onHP{i}.Qini;
            delMdeltau = onHP{i}.Qinidot;
            B = onHP{i}.Qini;
            delBdeltau = onHP{i}.Qinidot;
            
        else
            if i~=sys.Nhp-1
                
                if i == 1
                    M = recomposeSOSmatrix(Qvals(i,:), mM, sys.xdim-1); 
                    B = reshape(Bvals(i,:),sys.xdim-1,sys.xdim-1);
                    MN = M;
                    BN = B;
                else
                    M = Mp;
                    B = Bp;
                end
                
                Mp = recomposeSOSmatrix(Qvals(i+1,:), mM, sys.xdim-1);
                Bp = reshape(Bvals(i+1,:),sys.xdim-1,sys.xdim-1);

                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;
                delMdeltau = (Mp - M)/(taup - tau);
                delBdeltau = (Bp - B)/(taup - tau);

            else
                M = Mp;
                B = Bp;
                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;
                delMdeltau = (MN - M)/(taup - tau);
                delBdeltau = (BN - B)/(taup - tau);
            end
                
        end
        
        for j= 1:sys.uncer_count
            [sm1,cm1] = SOSmatrix(rho,degs.sm1_dU,degs.sm1_dL,sys.xdim-1);
            
            
            uncer_real = sys.permutvector(j,:) * sys.uncermax; 
            [Arho,taudotnum,rhodot] = computeTransverseJacobian(rho,sys,fns,onHP{i},uncer_real,numsets);
            
            Mrhodot = dotMetricMatrix(rho,M,rhodot,sys.xdim-1,numsets);
            C_cond = delMdeltau * taudotnum + Mrhodot + taudotden * (Arho' * M + M * Arho);
            neg_fac = clean(1e-4 * taudotden*eye(sys.xdim-1),numsets.clean_thresh);
            
            eq1a = y' * (- C_cond + neg_fac - sm1 * (1 - rho'* B *rho) ) * y;
            eq1b = y' * sm1 * y;
            
            F1 = [sos(eq1a),sos(eq1b)];
            [sol1] = solvesos(F1, [], numsets.sdpsetting1, cm1);
            prob = [prob,sol1.problem];
            if (prob==0)
                multi_mats{i}{j}.cm1 = value(cm1);
            end   
            
            [h1,hc1, ~] = polynomial(rho,degs.h1_dU,degs.h1_dL);
            
            Bdot = clean(rho'*(delBdeltau * taudotnum)*rho + 2*rho*B*rhodot,numsets.clean_thresh);
            
            eq4a = - Bdot + neg_fac + h1 * (1 - rho'* B *rho);
            
            F4 = [sos(eq4a)];
            [sol4] = solvesos(F4, [], numsets.sdpsetting1, hc1);
            prob = [prob,sol4.problem];
            if all(prob==0)
                multi_mats{i}{j}.hc1 = value(hc1)';
            end 
            
        end
        
        
        eq2 = taudotden  - s2 * (1 - rho'* B *rho);  
        eq3 = -((alphai - rho'*Smatfix*rho) * s3 + (rho'* B *rho - 1));
        
         
        F2 = [sos(eq2), sos(s2)]; 
        F3 = [sos(eq3), sos(s3)]; 
         
        [sol2] = solvesos(F2,[], numsets.sdpsetting1,c2);
        prob = [prob,sol2.problem];        
        [sol3] = solvesos(F3,[], numsets.sdpsetting1,c3);
        prob = [prob,sol3.problem];

        
        % in case feasibility test is infeasible 
        if any(prob~=0)
            fprintf(1,'In round %d the following errors were found:\n',i);
            disp(prob)
            error('Error in multiplier search')
        end
        
        % collect results to return
        multis{i}.c2 = value(c2)'; 
        multis{i}.c3 = value(c3)';

    end
    
    fprintf(1,'-------Multiplier search successfully finished-------\n\n')
end
                
    

function [alphares,Qvals,Bvals] = optimize_Metric(y,rho,fns,onHP,sys, numsets,degs,alpha_list,multis,multi_mats)

    Smatfix = numsets.ellipsemat;
    
    [~,~,mm1] = polynomial(rho,degs.sm1_dU,degs.sm1_dL);
    [~,~,m2] = polynomial(rho,degs.s2_dU,degs.s2_dL);
    [~,~,m3] = polynomial(rho,degs.s3_dU,degs.s3_dL);
    [~,~,mM] = polynomial(rho,degs.M_dU,degs.M_dL);
    [~,~,hm1] = polynomial(rho,degs.h1_dU,degs.h1_dL);

    F = []; 
    coeflist = [];%zeros(sys.Nhp-1,length(mM)*(sys.xdim-1)*(sys.xdim)/2);
    ballcoeflist = [];
    ObjFun  = [];
    
    for i=1:sys.Nhp-1
        
        if numsets.tvalpha == 1          
             alpha = sdpvar(1,1);
             alphatot = alpha_list(end,i) + alpha;
             alphaarray(i) = alpha;
             F = [F,alpha>=0];
             ObjFun = [ObjFun,-alpha];
        elseif numsets.tvalpha == 0
            if i == 1
                alpha = sdpvar(1,1);
                alphatot = alpha_list(end) + alpha;
                alphaarray = alpha;
                F = [F,alpha>=0];
                ObjFun = [ObjFun,-alpha];
            end
        end
       
        taudotden = onHP{i}.taudotden;
        if i < sys.Nhp-1           
            if i ==1
                [M,cM] = SOSmatrix([rho],degs.M_dU,degs.M_dL,sys.xdim-1);
                B = sdpvar(sys.xdim-1,sys.xdim-1);
                MN = M;
                BN = B;
                coeflist = [coeflist;cM];
                ballcoeflist = [ballcoeflist;reshape(B,1,(sys.xdim-1)^2)];
                F = [F,sos(y'*M*y)];
            else
                M = Mp;
                B = Bp;
            end

            [Mp,cMp] = SOSmatrix([rho],degs.M_dU,degs.M_dL,sys.xdim-1);
            Bp = sdpvar(sys.xdim-1,sys.xdim-1);
            coeflist = [coeflist;cMp];
            ballcoeflist = [ballcoeflist;reshape(Bp,1,(sys.xdim-1)^2)];
            F = [F,sos(y'*Mp*y)];
            
            tau = onHP{i}.tau;
            taup = onHP{i+1}.tau;
            delMdeltau = (Mp - M)/(taup - tau);  
            delBdeltau = (Bp - B)/(taup - tau); 
    
        else
            tau = onHP{i}.tau;
            taup = onHP{i+1}.tau;
            M = Mp;
            B = Bp;
            
            delMdeltau = (MN - M)/(taup - tau); 
            delBdeltau = (BN - B)/(taup - tau); 
        end

        for j = 1:sys.uncer_count
            sm1fix = recomposeSOSmatrix(multi_mats{i}{j}.cm1, mm1, sys.xdim-1);

            uncer_real = sys.permutvector(j,:) * sys.uncermax; 
            [Arho,taudotnum,rhodot] = computeTransverseJacobian(rho,sys,fns,onHP{i},uncer_real,numsets);
            
            Mrhodot = dotMetricMatrix(rho,M,rhodot,sys.xdim-1,numsets);
            C_cond = delMdeltau * taudotnum + Mrhodot + taudotden * (Arho' * M + M * Arho);
            neg_fac = clean(1e-4 * taudotden*eye(sys.xdim-1),numsets.clean_thresh);
            
            eq1 = y' * (- C_cond + neg_fac - sm1fix * (1 - rho'* B *rho) ) * y;
            
            F = [F,sos(eq1)];
            
            h1fix = multi_mats{i}{j}.hc1*hm1;
            
            Bdot = clean(rho'*(delBdeltau * taudotnum)*rho + 2*rho*B*rhodot,numsets.clean_thresh);
            eq4 = - Bdot + neg_fac + h1fix * (1 - rho'* B *rho);
            
            F = [F,sos(eq4)];
    
        end
        
        s2fix = clean(multis{i}.c2*m2,numsets.clean_thresh);
        s3fix = clean(multis{i}.c3*m3,numsets.clean_thresh);
        
        eq2 = taudotden  - s2fix * (1 - rho'* B *rho);  
        eq3 = -((alphatot - rho'*Smatfix*rho) * s3fix + (rho'* B *rho - 1));
        
        F = [F, sos(eq2), sos(eq3)];
    end

    [sol] = solvesos(F, ObjFun, numsets.sdpsetting2,[coeflist(:);ballcoeflist(:);alphaarray(:)]);
        
    if sol.problem ~= 0
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in Metric optimization')
    end

    Qvals = value(coeflist);
    Bvals = value(ballcoeflist);
    alphares = alpha_list(end,:) + value(alphaarray);
            
    fprintf('This is the alpha result: %f \n',max(alphares));

    fprintf(1,'-------Metric optimization successfully finished-------\n\n')

end


function [s,c] = SOSmatrix(rho,degMax,degMin,matrixsize)

    c = [];    
    for i = 1:matrixsize
        for j = i:matrixsize
            [sone,cone,~] = polynomial(rho,degMax,degMin);
            s(i,j) = [sone];
             if i~=j
                s(j,i) = [sone];
            end
            c =[c,cone'] ;
        end
    end   
          
     c = c; 
end

function [Mdot] = dotMetricMatrix(rho,M,fx,matrixsize,numsets)

       
    for i = 1:matrixsize
        for j = 1:matrixsize
            Mdot(i,j) = clean(jacobian(M(i,j),rho)*fx,numsets.clean_thresh);
        end
    end
       

end



function mat = recomposeSOSmatrix(crow, mcol, matrixsize)

    k = 0;
    len = length(mcol);
    for i = 1:matrixsize
        for j = i:matrixsize
            k = k+1;
            mat(i,j) = crow((k-1)*len+1:k*len)*mcol;
            if i~=j
                mat(j,i) = mat(i,j);
            end
        end
    end   
          
end


function [Arho,taudotnum,rhodot] = computeTransverseJacobian(rho,sys,fns,onHP,uncer_real,numsets)

    xx      = sdpvar(sys.xdim,1);
    Z       = onHP.Z;
    Zdot    = onHP.Zdot;
    v       = onHP.v;
    forb    = onHP.forb;
    xloc    = onHP.xloc;
    taudden = onHP.taudotden;
      
    fx_xx = fns.system_dyn_poly_uncertain(sys,onHP,xx,uncer_real);
    A = clean(jacobian(fx_xx,xx),numsets.clean_thresh);
    A = replace(A,xx,onHP.xloc);
    Arho = clean(Z*A*Z'+Zdot*Z'-Z*forb*(v'*A*Z'-v'*Zdot')/(v'*forb),numsets.clean_thresh);

    fx = fns.system_dyn_poly_uncertain(sys,onHP,xloc,uncer_real);
    taudotnum = clean(v'*fx,numsets.clean_thresh);
    
    rhodot = clean(Zdot*Z'*rho*taudotnum,numsets.clean_thresh) + clean(Z*fx*taudden,numsets.clean_thresh) - clean(Z*forb*taudotnum,numsets.clean_thresh); 
    rhodot = clean(rhodot,numsets.clean_thresh);

end





        

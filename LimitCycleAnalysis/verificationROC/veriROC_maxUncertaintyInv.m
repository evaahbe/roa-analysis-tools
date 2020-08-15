function veriROC_maxUncertaintyInv(sys,fns,filenames,pinit,rocnumsets,numsets)

    
    if ~exist(filenames.rocfinalresults,'file') 
        veriROC_maxROCinv(sys,fns,filenames,pinit,rocnumsets,[]);
    end
    
    numsets.ellipsemat = rocnumsets.ellipsemat;
    numsets.tvalpha = rocnumsets.tvalpha;
    rho = sdpvar(sys.xdim-1,1); % create n-1 dim state on hyperplane
    y = sdpvar(sys.xdim-1,1); % create n-1 dim indeterminant vector for matrix SOS check
    tau_array = linspace(0,sys.Tperiod,sys.Nhp)'; % create array of tau samples
    degs = numsets.degs;
    
    
    
    %initializing the variables
    try             % load results of already performed iteration steps
        load(filenames.uncerintermedresults)
    catch           % load the results from the ROC PRE maximization, then create intermediate results file for uncer for possible reinitialization
        load(filenames.rocfinalresults,'alpha','Qvals','gamma','Kvals','multis','multi_mats'); % the following loaded: alpha,Qvals,gamma,onHP.K,multis,multi_mats
        step    = 1;
        uncer_list(:,1) = -ones(sys.uncerdim,1);
        uncer_list(:,2) = sys.uncermax;
       
        try load(filenames.rocfinalresults,'Qvalsmet');
            Qvalsset = Qvals;
            Qvals = Qvalsmet;
            degs.B_dU = 0;
            degs.B_dL = 0;
        catch
            Qvalsset = Qvals;
            degs.B_dU = degs.M_dU;
            degs.B_dL = degs.M_dL;
        end
        save(filenames.uncerintermedresults,'step','Qvals','gamma','uncer_list','alpha','Kvals','Qvalsset','degs') 
    end
    
    
    
    numsets.alpha = alpha;
    
    for i = 1:sys.Nhp
        tau = tau_array(i);
             
        xorbi = ppval(sys.xorb, tau);
        [Zi, Zdoti] = projectOperator(sys, tau);
        xloci = xorbi + Zi' * rho; 
        vi = ppval(sys.v,tau);
        
        if sys.udim == 0 

            forbi = fns.system_dyn(xorbi,sys);  
            Kvals = [];
        else      
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
        onHP{i}.xloc        = xloci;
        onHP{i}.forb        = forbi;
        onHP{i}.xorb        = xorbi;
        onHP{i}.Z           = Zi;
        onHP{i}.Zdot        = Zdoti;
        onHP{i}.v           = vi;
        onHP{i}.Qvalsfix    = Qvalsset(i,:);
    end
         

    jj = length(size(uncer_list,2))-1;
    while jj<=numsets.iteration_max

        if (max(uncer_list(:,end) - uncer_list(:,end-1))) < numsets.convCrit*max(uncer_list(:,end-1))
            disp('Convergence reached!')
            break;
        end
        
        fprintf(1,'\n \n This is iteration nr.: %d \n',(size(uncer_list,2)-1))
            
        if step == 1   
            fprintf(1,'- step: Multiplier search\n\n')
            [multis, multi_mats] = find_multiplier(y,rho,fns,onHP,sys,numsets,degs,Qvals,uncer_list);
            step = 2;
            save(filenames.uncerintermedresults,'-append','step','multis','multi_mats')
        elseif step == 2 
            fprintf(1,'- step: Uncertainty bound maximization\n\n')
            [uncer_vec] = maximize_UncertaintyBounds(y,rho,fns, onHP, sys,uncer_list, degs,multi_mats,Qvals,numsets);
            uncer_list = [uncer_list,uncer_vec];
            if numsets.metric_search ==1
                step = 3;    
            else
                jj = jj+1;
                step = 1;
            end
            save(filenames.uncerintermedresults,'-append','uncer_list','step')
        elseif step == 3 % this step is only called if numsets.metric_search == 1
            fprintf(1,'- step: Metric search\n\n')
            [Qvals] = find_Metric(y,rho,fns, onHP, sys, degs, multis, multi_mats, uncer_list,numsets);
            step = 1;
            jj = jj+1;
            save(filenames.uncerintermedresults,'-append','Qvals','step')
        end
            
    end
    Qvals = Qvalsset;
    save(filenames.uncerfinalresults, 'multis','multi_mats','Qvals','step','alpha','gamma','Kvals','uncer_list','sys','numsets')  
end



function [multis,multi_mats] = find_multiplier(y,rho,fns,onHP,sys,numsets,degs,Qvals,uncer_list)


    [~,~,mM] = polynomial(rho,degs.M_dU,degs.M_dL);
    [~,~,mMfix] = polynomial(rho,degs.B_dU,degs.B_dL);
    uncer_vec = uncer_list(:,end);
    
    for i = 1:sys.Nhp-1
        prob = [];
        
        taudotden = onHP{i}.taudotden;    
        
        if numsets.tvalpha == 1          
            alphai = numsets.alpha(i);
        elseif numsets.tvalpha == 0      
            alphai = numsets.alpha;
        end
      
        if isempty(Qvals)
            error('Qvals should not be empty, take preliminary results from maxROC!')
            
        else
            if i~=sys.Nhp-1
                
                if i == 1
                    M = recomposeSOSmatrix(Qvals(i,:), mM, sys.xdim-1);
                    Mfix = recomposeSOSmatrix(onHP{i}.Qvalsfix, mMfix, sys.xdim-1);
                    MN = M;
                    MNfix = Mfix;
                else
                    M = Mp;
                end
                
                Mp = recomposeSOSmatrix(Qvals(i+1,:), mM, sys.xdim-1);
                Mpfix = recomposeSOSmatrix(onHP{i}.Qvalsfix, mMfix, sys.xdim-1);

                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;
                delMdeltau = (Mp - M)/(taup - tau);
                delMfixdeltau = (Mpfix - Mfix)/(taup - tau);
            else
                M = Mp;
                Mfix = Mpfix;
                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;
                delMdeltau = (MN - M)/(taup - tau);
                delMfixdeltau = (MNfix - Mfix)/(taup - tau);
            end
                
        end
        
        
        for j= 1:sys.uncer_count
            [sm1,cm1] = SOSmatrix(rho,degs.sm1_dU,degs.sm1_dL,sys.xdim-1);
            
            uncer_real = sys.permutvector(j,:) * uncer_vec; 
            [Arho,taudotnum,rhodot] = computeTransverseJacobian(rho,sys,fns,onHP{i},uncer_real,numsets);
            
            Mrhodot = dotMetricMatrix(rho,M,rhodot,sys.xdim-1,numsets);
            C_cond = delMdeltau * taudotnum + Mrhodot + taudotden * (Arho' * M + M * Arho);
            neg_fac = clean(1e-4 * taudotden*eye(sys.xdim-1),numsets.clean_thresh);
            
            
            eq1a = y' * (- C_cond +neg_fac - sm1 * (1 - rho'* Mfix *rho) ) * y;
            eq1b = y' * sm1 * y;
            
            F1 = [sos(eq1a),sos(eq1b)];
            [sol1] = solvesos(F1, [], numsets.sdpsetting1, cm1);
            prob = [prob,sol1.problem];
            if all(prob==0)
                multi_mats{i}{j}.cm1 = value(cm1);
            else
                fprintf(1,'In round %d the following errors were found:\n',i);
                disp(prob)
                error('Error in multiplier search')
            end
            
            [h1,hc1, ~] = polynomial(rho,degs.h1_dU,degs.h1_dL);
            Mfixrhodot = dotMetricMatrix(rho,Mfix,rhodot,sys.xdim-1,numsets);
            
            Vdot = clean(rho'*(delMfixdeltau * taudotnum + Mfixrhodot)*rho + 2*rho*Mfix*rhodot,numsets.clean_thresh);
            
            eq4a = - Vdot + neg_fac + h1 * (1 - rho'* Mfix *rho);
            
            F4 = [sos(eq4a)];
            [sol4] = solvesos(F4, [], numsets.sdpsetting1, hc1);
            prob = [prob,sol4.problem];
            if all(prob==0)
                multi_mats{i}{j}.hc1 = value(hc1)';
            else
                fprintf(1,'In round %d the following errors were found:\n',i);
                disp(prob)
                error('Error in multiplier search')
            end

        end
        
        
        multis = [];


    end
    
    fprintf(1,'-------Multiplier search successfully finished-------\n\n')
end
                
function [uncer_res] = maximize_UncertaintyBounds(y, rho,fns,onHP, sys, uncer_list, degs,multi_mats,Qvals,numsets)



    uncer_vec = sdpvar(sys.uncerdim,1);
    uncer_tot = uncer_list(:,end) + uncer_vec;
    
    [~,~,mm1] = polynomial(rho,degs.sm1_dU,degs.sm1_dL);
    [~,~,mM] = polynomial(rho,degs.M_dU,degs.M_dL);
    [~,~,mMfix] = polynomial(rho,degs.B_dU,degs.B_dL);
    [~,~,hm1] = polynomial(rho,degs.h1_dU,degs.h1_dL);
    
    F = [];
    for i=1:sys.Nhp-1
        if numsets.tvalpha == 1          
            alphai = numsets.alpha(i);
        elseif numsets.tvalpha == 0      
            alphai = numsets.alpha;
        end
        
        taudotden = onHP{i}.taudotden;
        
         if i~=sys.Nhp-1
                
                if i == 1
                    M = recomposeSOSmatrix(Qvals(i,:), mM, sys.xdim-1);
                    Mfix = recomposeSOSmatrix(onHP{i}.Qvalsfix, mMfix, sys.xdim-1);
                    MN = M;
                    MNfix = Mfix;
                else
                    M = Mp;
                    Mfix = Mpfix;
                end
                
                Mp = recomposeSOSmatrix(Qvals(i+1,:), mM, sys.xdim-1);
                Mpfix = recomposeSOSmatrix(onHP{i}.Qvalsfix, mMfix, sys.xdim-1);

                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;
                delMdeltau = (Mp - M)/(taup - tau);
                delMfixdeltau = (Mpfix - Mfix)/(taup - tau);
            else
                M = Mp;
                Mfix = Mpfix;
                tau= onHP{i}.tau;
                taup = onHP{i+1}.tau;
                delMdeltau = (MN - M)/(taup - tau);
                delMfixdeltau = (MNfix - Mfix)/(taup - tau);
            end
                
  
        for j = 1:sys.uncer_count
            sm1fix = recomposeSOSmatrix(multi_mats{i}{j}.cm1, mm1, sys.xdim-1);

            uncer_real = sys.permutvector(j,:) * uncer_tot; 
            [Arho,taudotnum,rhodot] = computeTransverseJacobian(rho,sys,fns,onHP{i},uncer_real,numsets);
            
            Mrhodot = dotMetricMatrix(rho,M,rhodot,sys.xdim-1,numsets);
            C_cond = delMdeltau * taudotnum + Mrhodot + taudotden * (Arho' * M + M * Arho);
            neg_fac = clean(1e-4 * taudotden*eye(sys.xdim-1),numsets.clean_thresh);
            
            eq1 = y' * (- C_cond + neg_fac - sm1fix * (1 - rho'* Mfix *rho) ) * y;
            
            F = [F,sos(eq1)];
            
            h1fix = clean(multi_mats{i}{j}.hc1*hm1,numsets.clean_thresh);
            Mfixrhodot = dotMetricMatrix(rho,Mfix,rhodot,sys.xdim-1,numsets);
            Vdot = clean(rho'*(delMfixdeltau * taudotnum + Mfixrhodot)*rho + 2*rho*Mfix*rhodot,numsets.clean_thresh);
            eq4 = - Vdot + neg_fac + h1fix * (1 - rho'* Mfix *rho);
            
            F = [F,sos(eq4)];
        end
             
    end
    
    F = [F, uncer_vec>=0];
    [sol] = solvesos(F, -uncer_vec, numsets.sdpsetting2,uncer_vec);
        
    if sol.problem ~= 0
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in maximization of uncertainty bounds.')
    end

    uncer_res = value(uncer_tot);
            
    fprintf('This is the uncertainty bounds result: %f \n',uncer_res);

    fprintf(1,'-------Maximization of uncertainty bounds successfully finished-------\n\n')
         
        


end


function [Qvals] = find_Metric(y, rho, fns,onHP, sys, degs, multis, multi_mats, uncer_list,numsets)

    Smatfix = numsets.ellipsemat;
    
    [~,~,mm1] = polynomial([rho],degs.sm1_dU,degs.sm1_dL);
    [~,~,mM] = polynomial(rho,degs.M_dU,degs.M_dL);
    [~,~,mMfix] = polynomial(rho,degs.B_dU,degs.B_dL);
    [~,~,hm1] = polynomial(rho,degs.h1_dU,degs.h1_dL);
    
    F = []; 
    coeflist = [];
    uncer_vec = uncer_list(:,end);

    for i=1:sys.Nhp-1
        
        if numsets.tvalpha == 1          
            alphai = numsets.alpha(i);
        elseif numsets.tvalpha == 0      
            alphai = numsets.alpha;
        end
       
        taudotden = onHP{i}.taudotden;
        if i < sys.Nhp-1           
            if i ==1
                [M,cM] = SOSmatrix([rho],degs.M_dU,degs.M_dL,sys.xdim-1);
                MN = M;
                coeflist = [coeflist;cM];
                F = [F,sos(y'*M*y)];
                Mfix = recomposeSOSmatrix(onHP{i}.Qvalsfix, mMfix, sys.xdim-1);
                MNfix = Mfix;
            else
                M = Mp;
                Mfix = Mpfix;
            end

            [Mp,cMp] = SOSmatrix([rho],degs.M_dU,degs.M_dL,sys.xdim-1);
            coeflist = [coeflist;cMp];
            F = [F,sos(y'*Mp*y)];
            
            tau = onHP{i}.tau;
            taup = onHP{i+1}.tau;
            delMdeltau = (Mp - M)/(taup - tau);  
            
            Mpfix = recomposeSOSmatrix(onHP{i}.Qvalsfix, mMfix, sys.xdim-1);
            delMfixdeltau = (Mpfix - Mfix)/(taup - tau);
        else
            tau = onHP{i}.tau;
            taup = onHP{i+1}.tau;
            M = Mp;
            Mfix = Mpfix;
            delMdeltau = (MN - M)/(taup - tau); 
            delMfixdeltau = (MNfix - Mfix)/(taup - tau);
        end
 
        
        for j = 1:sys.uncer_count
            sm1fix = recomposeSOSmatrix(multi_mats{i}{j}.cm1, mm1, sys.xdim-1);

            uncer_real = sys.permutvector(j,:) * uncer_vec; 
            [Arho,taudotnum,rhodot] = computeTransverseJacobian(rho,sys,fns,onHP{i},uncer_real,numsets);
            
            Mrhodot = dotMetricMatrix(rho,M,rhodot,sys.xdim-1,numsets);
            C_cond = delMdeltau * taudotnum + Mrhodot + taudotden * (Arho' * M + M * Arho);
            neg_fac = clean(1e-4 * taudotden*eye(sys.xdim-1),numsets.clean_thresh);
            
            eq1 = y' * (- C_cond +neg_fac - sm1fix * (1 - rho'* Mfix *rho) ) * y;
            
            F = [F,sos(eq1)];
            
            h1fix = clean(multi_mats{i}{j}.hc1*hm1,numsets.clean_thresh);
            Mfixrhodot = dotMetricMatrix(rho,Mfix,rhodot,sys.xdim-1,numsets);
            
            Vdot = clean(rho'*(delMfixdeltau * taudotnum + Mfixrhodot)*rho + 2*rho*Mfix*rhodot,numsets.clean_thresh);
            eq4 = - Vdot + neg_fac + h1fix * (1 - rho'* Mfix *rho);
            
            F = [F,sos(eq4)];

        end
    end

    [sol] = solvesos(F, [], numsets.sdpsetting3,[coeflist(:)]);
        
    if sol.problem ~= 0
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in Metric search')
    end

    Qvals = value(coeflist);
  
    fprintf(1,'-------Metric search successfully finished-------\n\n')

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





        

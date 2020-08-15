function recover_ROAstoch(sys,filenames,numsetsRE)


    load(filenames.finalresultsROAPCE,'Qvals','alpha','xiep')
    
    if sys.p ==0
        Q0 = Qvals;
        xiep0 = xiep;
        save(filenames.finalresultsROAstoch,'Q0','alpha','xiep0','sys','numsetsRE')
        return
    end

	degs = numsetsRE.degs;

    xi     = sdpvar((sys.p+1)*(sys.xdim),1);
    x0     = [xi(1)];
    xj     = [xi(2:sys.p+1)];
    xiep0  = [xiep(1)];
    for i = 1:sys.xdim-1
        x0    = [x0;xi((sys.p+1)*i+1)];
        xj    = [xj;xi((sys.p+1)*i+2:(sys.p+1)*(i+1))];
        xiep0 = [xiep0;xiep((sys.p+1)*i+1)];
    end

    mV = monolist([xi],degs.V_dU/2,degs.V_dL/2);
     
    if all(sys.varfix==0)
        
        mV0 = replace(mV,xj,zeros(size(xj)));
        mV0 = replace(mV0,x0,ones(size(x0)));
        mV0mat = mV0*mV0';
        Q0 = Qvals.*mV0mat;
        Q0 = Q0(:,any(Q0,1));
        Q0 = Q0(any(Q0,2),:);
        save(filenames.finalresultsROAstoch,'Q0','alpha','xiep0','sys','numsetsRE')
    else

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

        mx0 = monolist([x0],degs.V_dU/2,degs.V_dL/2);
        Q0 = diag(ones(length(mx0),1))/numsetsRE.initQ0scale;
        B0 = 10*diag(ones(length(x0),1))/numsetsRE.initQ0scale;
        B0matlist = [1e10,1e9];
        
        gamma = alpha;

        for jj=1:numsetsRE.iteration_max
            
            if B0matlist(end-1)- B0matlist(end) < numsetsRE.convCrit*B0matlist(end)
                break
            end

            fprintf(1,'\n \n This is iteration nr.: %d \n',jj)

            multis = find_multiplier(xi,x0,degs,Qvals,B0,Q0,sys,normsvec,alpha,gamma,numsetsRE);
            
            [Q0,B0] = optimize_Q0(xi,x0,degs,Qvals,multis,sys,normsvec,B0,alpha,gamma,numsetsRE);
            B0matlist = [B0matlist,det(B0)];

        end
        save(filenames.finalresultsROAstoch,'Q0','multis','B0','B0matlist','alpha','xiep0','sys','numsetsRE')
    end  
    
    
    
end     
 


function multis = find_multiplier(xi,x0,degs,Qvals,B0,Q0,sys,normsvec,alpha,gamma,numsetsRE)


    [s1,c1,~] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [s2,c2,~] = polynomial([x0],degs.s2_dU,degs.s2_dL);
    
    mx0 = monolist([x0],degs.V_dU/2,degs.V_dL/2);
    mV = monolist([xi],degs.V_dU/2,degs.V_dL/2);
   
    
    norms_mat = diag(normsvec);
    varfix_constraint = 0;
    i = 0;
    for k = 1:sys.xdim
        for l = k:sys.xdim
            if sys.varfix(k,l) ~= 0
                i = i+1;
                [hi,chi,~] = polynomial([xi],degs.hi_dU,degs.hi_dL);
                varfix_constraint = varfix_constraint + hi * (sys.varfix(k,l)-xi((sys.p+1)*(k-1)+1:(sys.p+1)*k)'*norms_mat*xi((sys.p+1)*(k-1)+1:(sys.p+1)*k));
                chi_array(i,:) = chi';
            end
        end
    end

    eq1 = - ((alpha-mx0'*Q0*mx0) * s1 - (alpha-mV'*Qvals*mV) + varfix_constraint);
    eq2 = - (gamma-x0'*B0*x0) * s2 + (alpha - mx0'*Q0*mx0);

    F1 = [sos(s1), sos(eq1)];
    F2 = [sos(eq2), sos(s2)];

    [sol1] = solvesos(F1, [], numsetsRE.sdpsetting1, [c1(:);chi_array(:)]);
    [sol2] = solvesos(F2, [], numsetsRE.sdpsetting1, c2);

    res1 = sol1.problem;
    res2 = sol2.problem;       
        
    if res1~=0 || res2~=0
        fprintf(1,'The following errors were found: %d %d \n',[res1,res2]);
        error('Error in multiplier search')
    else 
        multis.c1 = value(c1)'; 
        multis.c2 = value(c2)';
        multis.hi = value(chi_array);
    end
                            
    fprintf(1,'-------Multiplier search successfully finished-------\n\n')

end


function [Q0res,B0res] = optimize_Q0(xi,x0,degs,Qvals,multis,sys,normsvec,B0old,alpha,gamma,numsetsRE)
    

    [~,~,m1] = polynomial([xi],degs.s1_dU,degs.s1_dL);
    [~,~,m2] = polynomial([x0],degs.s2_dU,degs.s2_dL);
    [~,~,mh] = polynomial([xi],degs.hi_dU,degs.hi_dL);

    mx0 = monolist([x0],degs.V_dU/2,degs.V_dL/2);
    mV = monolist([xi],degs.V_dU/2,degs.V_dL/2);

    s1fix = multis.c1*m1;
    s2fix = multis.c2*m2;

    norms_mat = diag(normsvec);
    varfix_constraint = 0;
    i = 0;
    for k = 1:sys.xdim
        for l = k:sys.xdim
            if sys.varfix(k,l) ~= 0
                i = i+1;
                hifix = multis.hi(i,:)*mh;
                varfix_constraint = varfix_constraint + hifix * (sys.varfix(k,l)-xi((sys.p+1)*(k-1)+1:(sys.p+1)*k)'*norms_mat*xi((sys.p+1)*(k-1)+1:(sys.p+1)*k));
            end
        end
    end

    Q0 = sdpvar(length(mx0),length(mx0));
    B0diff = sdpvar(length(x0),length(x0));     
    B0 = B0old - B0diff;

    eq1 = -((alpha-mx0'*Q0*mx0) * s1fix - (alpha-mV'*Qvals*mV) + varfix_constraint);
    eq2 = -(gamma-x0'*B0*x0) * s2fix + (alpha-mx0'*Q0*mx0);

    F = [sos(eq1), sos(eq2), B0>=0, B0diff>=0, Q0>=0];
    [sol] = solvesos(F, -geomean(B0diff), numsetsRE.sdpsetting2, [Q0(:);B0diff(:)]);

    if sol.problem ~= 0
        fprintf(1,'The following error occured: %d\n',sol.problem);
        error('Error in Q0 optimization')
    end
    
    B0res = value(B0);
    Q0res = value(Q0);
    fprintf('This is the geomean result: %f \n',det(B0res));  
     
    fprintf(1,'-------Q0 optimization successfully finished-------\n\n')


end




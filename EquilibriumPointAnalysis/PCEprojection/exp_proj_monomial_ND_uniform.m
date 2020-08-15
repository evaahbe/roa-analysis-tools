function mono = exp_proj_monomial_ND_uniform(var, vari, fix_var, fix_vari,  monobase, I, obj, norm_dec)

    %the i-th row of vari contains the coefficients of the i-th entry of
    %the var vector.
    %norm_dec is 1 when coefficient dynamics are computed. It is zero when
    %subterms are calculated which do not get projected onto the basis
    %polys directly (e.g. as is the case in the taudot comutations of LC)
    %fix_vari are for example the mu_coefs
    
  
    if value(monobase) == 1
        if I==1 % constant variables only appear in the zeroth expansion coefficient equations!
            mono = 1;
        else
            mono = 0;
        end
        return
    end
    
    %for the vars (x,u..)
    j = 0;
    powervec = [];
    for i = 1:length(var)
        varexpnr = exponents(monobase,var(i));
        if varexpnr ~= 0
            j = j+1;
            actvar(j,:) = vari(i,:); %select only the variables actually appearing in the monomial, as PCE coefficients
            powervec(j) = varexpnr; % save the exponent of the actually appearing variables
        end
    end
    
    %for the fix_vars (mu,..)
    j = 0;
    fixpowervec = [];
    for i = 1:length(fix_var)
        fixvarexpnr = exponents(monobase,fix_var(i));
        if fixvarexpnr ~= 0
            j = j+1;
            fixactvar(j,:) = fix_vari(i,:); %select only the variables actually appearing in the monomial, as PCE coefficients
            fixpowervec(j) = fixvarexpnr; % save the exponent of the actually appearing variables
        end
    end
    
    dim = sum(powervec) + sum(fixpowervec);
    
    if ~isempty(powervec)
        %initialize the coefficient matrix 
        coef_mat = actvar(1,:); %the first one has to a row vector to build the kronecker product later!
        powervec(1) = powervec(1)-1;
    else
        coef_mat = fixactvar(1,:); %the first one has to a row vector to build the kronecker product later!
        fixpowervec(1) = fixpowervec(1)-1;
    end

    for i = 1:length(powervec)
        for j = 1:powervec(i)
            coef_mat = kron(actvar(i,:)',coef_mat);
        end
    end
        
    for i = 1:length(fixpowervec)
        for j = 1:fixpowervec(i)
            coef_mat = kron(fixactvar(i,:)',coef_mat);
        end
    end
        
    
    if norm_dec == 1
        basis_tensor = create_tensor_reshaped_Legendre((I-1),obj.p,dim);
    elseif norm_dec == 0
        basis_tensor = create_tensor_reshaped_Legendre_nonorm((I-1),obj.p,dim);
    else
        error('decide on normalization of tensor!')
    end

    mono = sum(sum(coef_mat.*basis_tensor));
        
    
end
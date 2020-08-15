function fxi = castPCEdynamics(vartot, varitot,sys,fxnom)

    norm_dec =1;
    var = [];
    vari = [];
    fix_var = [];
    fix_vari = [];

    for i = 1:length(varitot(:,1))
        if isnan(double(varitot(i,1)))
            var = [var;vartot(i)];
            vari = [vari;varitot(i,:)];
        else
            fix_var = [fix_var;vartot(i)];
            fix_vari = [fix_vari;varitot(i,:)];
        end
    end

    
    if strcmp(sys.distType,'uniform')

        for k = 1:length(fxnom)
            [coefs,monobase] = coefficients(fxnom(k),vartot);

            for i =1:sys.p+1 %projection

                fxicomp = 0;
                for j = 1:length(monobase)
                    mono = exp_proj_monomial_ND_uniform(var,vari,fix_var, fix_vari,monobase(j),i,sys,norm_dec);
                    fxicomp = fxicomp+coefs(j)*mono;
                end

                indi = (k-1)*(sys.p+1)+i;
                fxi(indi) = fxicomp;

            end
        end
        
    elseif strcmp(sys.distType,'normal')
    
        for k = 1:length(fxnom)
            [coefs,monobase] = coefficients(fxnom(k),vartot);

            for i =1:sys.p+1 %projection

                fxicomp = 0;
                for j = 1:length(monobase)
                    mono = exp_proj_monomial_ND_normal(var,vari,fix_var, fix_vari,monobase(j),i,sys,norm_dec);
                    fxicomp = fxicomp+coefs(j)*mono;
                end

                indi = (k-1)*(sys.p+1)+i;
                fxi(indi) = fxicomp;

            end
        end
    else
        error('Distribution Type unknown.')
    end


end
function [eta,vispline,vidotspline] = hyperplane_base_auto(sys,fns)
    
    tvec = sys.tvec;
    epsi = 1e-5;
    dt = 0.000001/2;
    
    fx      = zeros(sys.xdim,length(tvec));
    vclass  = zeros(sys.xdim,length(tvec));
    xloc    = zeros(sys.xdim,length(tvec));
    
    %classical v-vectors
    for k = 1:length(tvec)
        xloc(:,k) = ppval(sys.xorb,tvec(k));
        fx(:,k) = fns.system_dyn(xloc(:,k),sys);
        vclass(:,k) = fx(:,k)./sqrt(sum(fx(:,k).^2));
    end
        
        
    if sys.trafo == 1
        
        if sys.xdim == 2
            eta = 0;
        else
            %construct coordinate system of hyperplane
            eta_one = mean(vclass,2)/norm(mean(vclass,2));

            for k = 1:length(tvec)
                if (sum(vclass(:,k).*eta_one) <=(1+epsi)) && (sum(vclass(:,k).*eta_one) >=(1-epsi))
                    error('numerically problematic eta base vector')
                end
            end

            eta_is = null(eta_one');
            eta = [eta_one, eta_is];

            if rank(eta)~=sys.xdim
                error('Something is wrong with the hyperplane base vectors')
            end
        end
        
        % construct v-vector spline
        vispline = spline(tvec,vclass);  
           
        % construct dvdtau-vector spline  
        vidot = zeros(sys.xdim,length(tvec));
        if isfield(fns,'system_dyn_deriv')  % --> with dynamics derivative, if available 
            for k = 1:length(tvec)
                dfxdt = fns.system_dyn_deriv(xloc(:,k),sys);
                vidot(:,k) = dfxdt./sqrt(sum(fx(:,k).^2)) + fx(:,k)./(-(sum(fx(:,k).^2))^(3/2)) .* dot(fx(:,k),dfxdt);
            end
        else                            % --> numerically, otherwise
            for k = 1:length(tvec)
                xa = ppval(sys.xorb,tvec(k)+dt);
                xb = ppval(sys.xorb,tvec(k)-dt); 
                fxa = fns.system_dyn(xa,sys);
                fxb = fns.system_dyn(xb,sys);
                va = fxa./sqrt(sum(fxa.^2));
                vb = fxb./sqrt(sum(fxb.^2));
                vidot(:,k) = (va-vb)./(2*dt);
            end
        end
  
        vidotspline = spline(tvec,vidot);
    
    elseif sys.trafo == 2
        
        vc = zeros(sys.xdim,length(tvec));
        %construct vectors from trajectory point to center point
        for k = 1:length(tvec)
            xtraj = ppval(sys.xorb,tvec(k));
            vc(:,k) = -(-xtraj+sys.xc)/(sqrt(sum((-xtraj+sys.xc).^2)));
            if (sum(vclass(:,k).*vc(:,k)) <=(1+epsi)) && (sum(vclass(:,k).*vc(:,k)) >=(1-epsi))
                error('numerically problematic or unsuitable center point!')
            end
%             vi = -[-vc(2,k);vc(1,k)];
%             if sum(vclass(:,k).*vi)<=0
%                 disp(k)
%                 stop = 1;
%             end
        end

        if sys.xdim == 2
            vi = [-vc(2,:);vc(1,:)]; %vc is projection operator here for xdim=2!
            if sum(vclass(:,1).*vi(:,1))<=0 
                vi = -[-vc(2,:);vc(1,:)];
            end     
            eta = 0;
        elseif sys.xdim >=3
            rotv = zeros(sys.xdim,length(tvec));
            for k = 1:length(tvec)
                %rotv(:,k) = null([vclass(:,k),vc(:,k)]');
                rotv(:,k) = cross(vclass(:,k),vc(:,k));
            end 
            fixv = mean(rotv,2)/norm(mean(rotv,2));
            %do the check
            for k = 1:length(tvec)
                if (sum(vclass(:,k).*fixv) <=(1+epsi)) && (sum(vclass(:,k).*fixv) >=(1-epsi))
                    error('numerically problematic or unsuitable fixed vector!')
                end
                if (sum(vc(:,k).*fixv) <=(1+epsi)) && (sum(vc(:,k).*fixv) >=(1-epsi))
                    error('numerically problematic or unsuitable fixed vector!')
                end
            end
            %compute the desired v vector
            vvec = zeros(sys.xdim,length(tvec));
            vi = zeros(sys.xdim,length(tvec));
            for k = 1:length(tvec)
                %vvec(:,k) = null([vc(:,k),fixv]');
                vvec(:,k) = cross(vc(:,k),fixv);
                vi(:,k) = vvec(:,k)./norm(vvec(:,k));
            end
            
            eta_one = mean(vvec,2)/norm(mean(vvec,2));
            eta_is = null(eta_one');
            eta = [eta_one, eta_is];

            if rank(eta)~=sys.xdim
                error('Something is wrong with the hyperplane base vectors')
            end
        
        end
     
         % construct v-vector spline
        vispline = spline(tvec,vi);
              
        % construct dvdtau-vector spline (numerically)
        vidot = zeros(sys.xdim,length(tvec));
        for k = 1:length(tvec)    
            va = ppval(vispline,tvec(k)-dt);
            vb = ppval(vispline,tvec(k)+dt);
            vidot(:,k) = (vb-va)./(2*dt);
        end
        
        vidotspline = spline(tvec,vidot);
        
    else 
        error('Transverse coordinate transformation not specified')
    end

   
end
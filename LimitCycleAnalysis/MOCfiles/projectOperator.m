
% function computes projection operator as in eq.VI.1.2 and VI.1.3, p. 215/216, Ordinary Differential Equations, Hale (1980) 

function [Zi, Zidot] = projectOperator(sys,tau)
 
        vi_tau = ppval(sys.v,tau);
        vidot_tau = ppval(sys.vdot,tau);

        if sys.xdim == 2
            Zi = [-vi_tau(2), vi_tau(1)];
            Zidot = [-vidot_tau(2), vidot_tau(1)];
        else
            eta = sys.eta;

            Zi    = zeros(sys.xdim,sys.xdim-1);
            Zidot = zeros(sys.xdim,sys.xdim-1);

            for j = 2:sys.xdim
                Zi(:,j-1) = eta(:,j)-(eta(:,j)'*vi_tau)/(1+eta(:,1)'*vi_tau).*(eta(:,1)+vi_tau);
                
                Zidot(:,j-1) = -(((eta(:,j)' * vidot_tau) .* (eta(:,1) + vi_tau)...
                    + (eta(:,j)' * vi_tau) .* vidot_tau) .* (1+eta(:,1)' * vi_tau) ...
                    -(eta(:,j)' *vi_tau) .* (eta(:,1) + vi_tau) .* (eta(:,1)'* vidot_tau))/(1 + eta(:,1)' * vi_tau)^2;
            end

            Zi = Zi';
            Zidot = Zidot';
        end
        


end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
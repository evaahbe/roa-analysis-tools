    
function [epsimu] = findPCEEPsimu(xi,fxi,sys)
    
    epfinding_threshold = 1e-4;

    syms xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)
    if ((sys.p+1)*sys.xdim) == 1
        xs = [xs1];
        fxi_sym = makesdpvarsym(fxi,xi,xs); 
        fxi_sym_t = subs(fxi_sym,[xs1],[xst1(t)]);
        xst = [xst1(t)];
    elseif ((sys.p+1)*sys.xdim) == 2
        xs = [xs1;xs2];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2],[xst1(t) xst2(t)]);
        xst = [xst1(t) xst2(t)];
    elseif ((sys.p+1)*sys.xdim) == 3
        xs = [xs1;xs2;xs3];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3],[xst1(t) xst2(t) xst3(t)]);
        xst = [xst1(t) xst2(t) xst3(t)];
    elseif ((sys.p+1)*sys.xdim) == 4
        xs = [xs1;xs2;xs3;xs4];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4],[xst1(t) xst2(t) xst3(t) xst4(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t)];
    elseif ((sys.p+1)*sys.xdim) == 5
        xs = [xs1;xs2;xs3;xs4;xs5];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t)];
    elseif ((sys.p+1)*sys.xdim) == 6
        xs = [xs1;xs2;xs3;xs4;xs5;xs6];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t)];
    elseif ((sys.p+1)*sys.xdim) == 7
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t)];
    elseif ((sys.p+1)*sys.xdim) == 8
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7;xs8];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t)];
    elseif ((sys.p+1)*sys.xdim) == 9
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7;xs8;xs9];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t)];
    elseif ((sys.p+1)*sys.xdim) == 10
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7;xs8;xs9;xs10];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)];
    else
        error('System and PCE dimensions not supported')
    end
    
    fxi_sym_fun = odeFunction(fxi_sym_t,xst);
    
    xi0 = kron(sys.x0,[1;zeros(sys.p,1)]);
    tinterval = linspace(0,sys.Tintegral,20000);
    [~, xvec] = ode45(fxi_sym_fun, tinterval,xi0);
        
    ep_cand = xvec(end,:);
   
    if sum(abs(xvec(end-10:end,:)-ep_cand))<epfinding_threshold
        epsimu = ep_cand;
    else
        fprintf(1,'Consider the following warning:')
        warning('Ran into problems with PCE EP simulation, trying larger integration time..')
        T = 1000;
        tinterval = linspace(0,T,50000);
        [tvec, xvec] = ode45(fxi_sym_fun, tinterval,xi0);
        ep_cand = xvec(end,:);
        if sum(abs(xvec(end-10:end,:)-ep_cand))<epfinding_threshold
            epsimu = ep_cand;
        else
            figure
            hold on
            for i = 1:sys.xdim*(sys.p+1)
                plot(tvec,xvec(:,i))
            end
            error('PCE EP could not be found with given integration settings (see figure)! Try changing: 1. initial condition (sys.x0), ode solver type (default:ode45), accuracy of PCE EP finding threshold')
         end
   end
end


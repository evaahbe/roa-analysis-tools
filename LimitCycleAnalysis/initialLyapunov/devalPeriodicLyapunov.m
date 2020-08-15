function [P, Pdot] = devalPeriodicLyapunov(tvec,sys,fns,pinit)
    
    Q = pinit.Qmat;
    Q0 = pinit.Q0mat;

    rounds = sys.rounds;
    
    tint = [tvec];
    
    for i = 2:1:rounds 
        tint = [tint;tvec(2:end)+(i-1)*tvec(end)]; 
    end
    
    for i = 1:length(tint)
        tt = tint(i);
        [Zi(:,:,i),Zdoti(:,:,i)] = projectOperator(sys,tt);
    end
        
    sys.zispline = spline(tint,Zi);
    sys.zidotspline = spline(tint,Zdoti);
    
    tint = fliplr(tint'); %integration starts at end, so flip time 
    %tintext = linspace(tint(1),tint(end),20000);
    
    
    
    [Z,~] = projectOperator(sys,tint(1));
    Q0 = Z*Q0*Z';
    
    
    for j = 1:length(tvec)
        tau = tvec(j);   
        P{j} = deval_matrix_ode(@ode45,@(t,P) periodicLyapunov( t, P, Q, sys, fns), tint, Q0, tau); 
        Pdot{j} = periodicLyapunov(  tau, P{j}, Q, sys, fns);
    end

end


function Pdot = periodicLyapunov(t, P, Q, sys, fns)
    %[Z,Zdot] = projectOperator(sys,t);
    Z = fnval(sys.zispline,t);
    Zdot = fnval(sys.zidotspline,t);
    A = fns.system_linearization(sys,t);
    v = ppval(sys.v,t);
    xorb = ppval(sys.xorb,t);
    forb = fns.system_dyn(xorb,sys);
    Qtran = Z*Q*Z';
    Atran = Z * A * Z' + Zdot*Z' - Z*forb*(v'*A*Z'-v'*Zdot')/(v'*forb);
    Pdot = -(Atran.'*P + P*Atran + Qtran); 
end
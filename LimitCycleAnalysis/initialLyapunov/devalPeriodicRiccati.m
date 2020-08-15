function [P, Pdot,K] = devalPeriodicRiccati(tvec,sys,fns,pinit)
    
    Q = pinit.Qmat;
    Q0 = pinit.Q0mat;
    R = pinit.R;
    
    rounds = sys.rounds;
    
    tint = [tvec];
    

    for i = 2:1:rounds 
        tint = [tint;tvec(2:end)+(i-1)*tvec(end)]; 
    end
    
    tint = fliplr(tint'); %integration starts at end, so flip time
    
    [Z,~] = projectOperator(sys,tint(1));
    Q0 = Z*Q0*Z';
    
    for j = 1:length(tvec)
        tau = tvec(j);
        
        P{j} = deval_matrix_ode(@ode45,@(t,P) periodicRiccati( t, P, Q, R, sys, fns), tint, Q0, tau); 
        
        Pdot{j} = periodicRiccati(  tau, P{j}, Q, R, sys, fns);
        
        [Z,~] = projectOperator(sys,tau);
        [~,B] = fns.system_linearization(sys,tau);
        v = ppval(sys.v,tau);
        xorb = ppval(sys.xorb,tau);
        uorb = ppval(sys.uorb,tau);
        forb = fns.system_dyn(xorb,uorb,sys);
        Btran = Z * B - Z*forb*v'*B/(v'*forb);
        K{j} = -R^(-1) *Btran.'* P{j};             
    end

end


function Pdot = periodicRiccati(t, P, Q, R, sys, fns)
    [Z,Zdot] = projectOperator(sys,t);
    [A,B] = fns.system_linearization(sys,t);
    v = ppval(sys.v,t);
    xorb = ppval(sys.xorb,t);
    uorb = ppval(sys.uorb,t);
    forb = fns.system_dyn(xorb,uorb,sys);
    Qtran = Z*Q*Z';
    Atran = Z * A * Z' + Zdot*Z' - Z*forb*(v'*A*Z'-v'*Zdot')/(v'*forb);
    Btran = Z * B - Z*forb*(v'*B/(v'*forb));
    Pdot = -(Atran.'*P + P*Atran - P*Btran*R^(-1)*Btran.'*P + Qtran); 
end
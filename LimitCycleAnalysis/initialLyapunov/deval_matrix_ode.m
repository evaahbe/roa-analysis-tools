function solattau = deval_matrix_ode(solver,functiontosolve,tint,Qf,tau)


    % reshape into a vector
    x0 = reshape(Qf,size(Qf,1)*size(Qf,2),1);
    solattau= deval(solver(@(t,x) reshaper_fun(t,x,functiontosolve,Qf),tint,x0),tau);
    solattau = reshape(solattau,size(Qf));



end
  
function vecdot = reshaper_fun(t,vec,functiontosolve,Qf)
    mat = reshape(vec,size(Qf));
    matdot = functiontosolve(t,mat);
    vecdot = reshape(matdot,size(matdot,1)*size(matdot,2),1);
end
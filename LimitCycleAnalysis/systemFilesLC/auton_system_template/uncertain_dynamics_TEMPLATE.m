function fx = uncertain_dynamics_TEMPLATE(sys,onHP,xloc,uncer_real)
    
    x(1) = xloc(1);
    x(2) = xloc(2);
    mu = sys.mu;
    
    % fill in the appropriate uncertain dynamics in poly approx if not poly
    % already
    % insert (mu(i)+uncer_real(i)) for every uncertain parameter!
    x1dot = 1;
    x2dot = (mu+uncer_real) * 1;
    
    fx = [x1dot;x2dot];
end
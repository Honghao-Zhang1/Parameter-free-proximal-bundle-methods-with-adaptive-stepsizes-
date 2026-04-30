function  [f,g]=helper_dual_oracle(q,W,T,h,pi,x)
    % names = { 'y'; 'u'};
    model.A = sparse([W,T]);
    model.obj = [q; -pi];
    model.rhs = h ;
    % model.lb= 0 ;
    % model.ub= 1 ;
    model.sense = '>';
    % type = 
    model.vtype = 'B';
    % model.varnames = names;

    % gurobi_write(model, 'mip1.lp');

    params.outputflag = 0;

    result = gurobi(model, params);
    % y = result.x(1:120);      
    u = result.x(121:240);
    g = -u+x;
    % pi'*x
    f = result.objval+pi'*x;
end
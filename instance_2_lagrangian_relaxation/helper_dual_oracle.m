function  [f,g]=helper_dual_oracle(q,W,T,h,pi,x)
    q  = q(:);
    pi = pi(:);
    x  = x(:);

    n = numel(pi);

    model.A = sparse([W,T]);
    model.obj = [q; -pi];
    model.rhs = h;
    model.sense = '>';

    % y is binary, u is continuous in [0,1]^n
    model.lb = zeros(2*n,1);
    model.ub = ones(2*n,1);
    model.vtype = [repmat('B', n, 1); repmat('C', n, 1)];

    params.outputflag = 0;

    result = gurobi(model, params);

    u = result.x(n+1:2*n);

    g = -u + x;
    f = result.objval + pi' * x;
end

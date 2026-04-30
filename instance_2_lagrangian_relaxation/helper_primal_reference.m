  function opt = helper_primal_reference(q,W,T,h,x)
    % names = { 'y'};

    model.A = sparse(W);
    model.obj = q;
    model.rhs = h - T*x;
    model.sense = '>';
    model.vtype = 'B';
    % model.varnames = names;

    % gurobi_write(model, 'mip1.lp');

    params.outputflag = 0;

    result = gurobi(model, params);
     opt = result.objval;
     % result.Pi
    % if strcmp(result.status, 'OPTIMAL');
    %   fprintf('Optimal objective: %e\n', result.objval);
    %   disp(result.x)
    % else
    %   fprintf('Optimization returned status: %s\n', result.status);
    % end
end
function x = helper_original_problem(c,d,A,C,b,Bundle_pi,Bundle_L,k)
    model.A = zeros(k+50,241);
    model.rhs = zeros(k+50,1);
    for i = 1:k
        % pi = Bundle_pi(i,:);
        model.A(i,:) = [ -Bundle_pi(i,:) zeros(1,120) 1];
        model.rhs(i,1) = Bundle_L(i);
    end
    model.A(k+1:k+50,:)= sparse([A,C,zeros(50,1)]);
    model.A = sparse(model.A);
    model.obj = [c', d' ,1];
    model.rhs(k+1:k+50,1)= b';
    % size(model.rhs)
    model.sense = '>';
    % type = ;
    % for i=1:240
    %     type(i) = 'B';
    % end
    % type(241)='C';
    % type
    model.vtype = char(ones(241,1) * 'B');
    model.vtype(end) = 'C';
    % model.varnames = names;

    % gurobi_write(model, 'mip1.lp');

    params.outputflag = 0;

    result = gurobi(model, params);
    x = result.x(1:120);
end 
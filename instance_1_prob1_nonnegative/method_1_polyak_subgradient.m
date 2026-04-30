function stats = method_1_polyak_subgradient(x0, eps_target, A, b, time_limit)
%METHOD_1_POLYAK_SUBGRADIENT Pol-Sub.
%
%   stats = method_1_polyak_subgradient(x0, eps_target, A, b, time_limit)
%   runs the projected subgradient method with Polyak stepsize on
%
%       min ||A*x-b||_1  subject to x >= 0.
%
%   This is the Pol-Sub baseline for the nonnegative l1 feasibility
%   experiment. It uses the Polyak stepsize lambda = f(x)/||g(x)||^2; for
%   these generated feasibility instances, phi* = 0.

    x = x0;
    iter = 0;
    t_start = tic;
    fbest = norm(A*x' - b, 1);

    while fbest > eps_target
        g = A' * sign(A*x' - b);
        fx = norm(A*x' - b, 1);
        lam = fx / max(norm(g)^2, 1e-16);
        x = max(x - lam*g', 0);

        if fx < fbest
            fbest = fx;
        end

        iter = iter + 1;
        if toc(t_start) > time_limit
            fprintf('ratio = %.4f\n', fbest / max(eps_target, 1e-16));
            disp('----out of time----')
            break
        end
    end

    stats = struct( ...
        'x', x, ...
        'fbest', fbest, ...
        'iterations', iter, ...
        'cpu_time', toc(t_start), ...
        'hit_time_limit', toc(t_start) > time_limit);
end

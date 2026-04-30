function [x, stats] = method_7_apl(x, eps, A, b, S, opts)
%METHOD_7_APL Accelerated proximal level comparison baseline.
%
%   [x, stats] = method_7_apl(x, eps, A, b, S) approximately solves
%
%       min ||A*x - b||_1  subject to x >= 0.
%
%   S is the maximum number of outer APL phases. Optional parameters can be
%   passed through opts:
%
%       opts.beta       Level parameter in (0,1). Default: 0.5.
%       opts.theta      Phase-reduction parameter in (0,1). Default: 0.5.
%       opts.max_inner  Maximum inner iterations per phase.
%       opts.time_limit Time limit in seconds.
%
%   stats records the exit reason, CPU time, phase count, and whether any
%   phase reached its inner iteration cap.
%
%   APL is kept as an external accelerated proximal-level comparison
%   baseline for the same nonnegative l1 feasibility model. Its projection
%   step uses helper_project_nonnegative_two_cuts.m from this folder.

    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    opts = fill_l1_nonnegative_apl_options(opts);

    t_start = tic;
    stats = init_l1_nonnegative_apl_stats(opts.max_inner);

    x = max(x(:), 0);
    b = b(:);

    hit_inner_cap = false;
    for s = 1:S
        [x, phase_info] = APL_reduction_prob1( ...
            x, eps, A, b, opts.beta, opts.theta, opts.max_inner, t_start, opts.time_limit);
        hit_inner_cap = hit_inner_cap || phase_info.hit_cap;

        stats.phases_run = s;
        stats.last_phase_iterations = phase_info.iterations;
        stats.last_phase_improve = phase_info.improve;
        stats.any_hit_inner_cap = hit_inner_cap;
        stats.num_phases_hit_inner_cap = stats.num_phases_hit_inner_cap + double(phase_info.hit_cap);

        Ub = norm(A * x - b, 1);
        if phase_info.hit_time_limit
            stats.exit_reason = 'time_limit';
            stats.cpu_time = toc(t_start);
            fprintf('ratio = %.4f\n', Ub / max(eps, 1e-16));
            disp('----out of time----');
            return
        end

        if Ub < eps
            stats.exit_reason = 'eps';
            stats.cpu_time = toc(t_start);
            fprintf('Reached eps: s=%d, Ub=%.6e, cpu=%.3f\n', s, Ub, stats.cpu_time);
            return
        end

        if toc(t_start) > opts.time_limit
            stats.exit_reason = 'time_limit';
            stats.cpu_time = toc(t_start);
            fprintf('ratio = %.4f\n', Ub / max(eps, 1e-16));
            disp('----out of time----');
            return
        end
    end

    stats.exit_reason = 'completed_S';
    stats.cpu_time = toc(t_start);
end

function [xupper, info] = APL_reduction_prob1(x, eps, A, b, beta, theta, max_inner, t_start, time_limit)
    k = 1;

    x = max(x(:), 0);
    center = x;
    xupper = x;

    Ub = norm(A * xupper - b, 1);
    Ub0 = Ub;

    % LB = 0, so l = (1-beta)*Ub0.
    l = (1 - beta) * Ub0;
    threshold = (1 - theta) * l + theta * Ub0;

    while Ub > threshold
        if toc(t_start) > time_limit
            info.iterations = k - 1;
            info.hit_cap = false;
            info.hit_time_limit = true;
            info.improve = max(0, Ub0 - Ub);
            return
        end

        alpha = 2 / (k + 1);
        xlower = (1 - alpha) * xupper + alpha * x;

        r = A * xlower - b;
        g = A' * sign(r);
        fvalue = norm(r, 1) - g' * xlower;
        rhs = l - fvalue;

        a = x - center;
        if norm(a) < 1e-12
            x = proj_Rplus_onecut(center, g, rhs, 'le', 1e-8);
        else
            x = helper_project_nonnegative_two_cuts(center, g, rhs, a, a' * x, 1e-8);
        end

        xtemp = (1 - alpha) * xupper + alpha * x;
        ftemp = norm(A * xtemp - b, 1);

        if ftemp < Ub
            xupper = xtemp;
            Ub = ftemp;
        end

        if Ub < eps
            fprintf('Inner reached eps: k=%d, Ub=%.6e, cpu=%.3f\n', k, Ub, toc(t_start));
            break
        end

        k = k + 1;
        if k > max_inner
            break
        end
    end

    info.iterations = min(k, max_inner);
    info.hit_cap = (k > max_inner) && (Ub > threshold);
    info.hit_time_limit = false;
    info.improve = max(0, Ub0 - Ub);
end

function stats = init_l1_nonnegative_apl_stats(max_inner)
    stats = struct( ...
        'exit_reason', '', ...
        'cpu_time', 0, ...
        'phases_run', 0, ...
        'last_phase_iterations', 0, ...
        'last_phase_improve', 0, ...
        'any_hit_inner_cap', false, ...
        'num_phases_hit_inner_cap', 0, ...
        'max_inner', max_inner);
end

function opts = fill_l1_nonnegative_apl_options(opts)
    defaults.beta = 0.9;
    defaults.theta = 0.35;
    defaults.max_inner = 200;
    defaults.time_limit = 3600;

    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end
end

function x = proj_Rplus_onecut(x0, c, d, sense, tol)
    x0 = x0(:);
    c = c(:);
    proj = @(z) max(z, 0);

    x = proj(x0);
    val = c' * x;

    if strcmpi(sense, 'le')
        if val <= d + tol
            return
        end
        lo = 0;
        hi = 1;
        while c' * proj(x0 - hi * c) > d + tol
            hi = 2 * hi;
        end
        for k = 1:60
            mid = (lo + hi) / 2;
            xm = proj(x0 - mid * c);
            if c' * xm <= d + tol
                hi = mid;
                x = xm;
            else
                lo = mid;
            end
        end
    else
        if val >= d - tol
            return
        end
        lo = 0;
        hi = 1;
        while c' * proj(x0 + hi * c) < d - tol
            hi = 2 * hi;
        end
        for k = 1:60
            mid = (lo + hi) / 2;
            xm = proj(x0 + mid * c);
            if c' * xm >= d - tol
                hi = mid;
                x = xm;
            else
                lo = mid;
            end
        end
    end
end

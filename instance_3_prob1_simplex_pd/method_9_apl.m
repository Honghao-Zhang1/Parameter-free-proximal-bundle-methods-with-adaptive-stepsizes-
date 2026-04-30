function [xbest, fbest, LB, stats] = method_9_apl(x, eps, A, b, S, D, max_inner, opts)
    % Accelerated proximal level method for
    %   min ||A*x - b||_1  s.t. x >= 0, sum(x) <= D
    %
    % Uses:
    % - one fresh level cut
    % - one Euclidean prox halfspace
    % - lower-bound updates from affine minorants over the bounded set

    if nargin < 7 || isempty(max_inner)
        max_inner = [];
    end
    if nargin < 8 || isempty(opts)
        opts = struct();
    end

    if isstruct(max_inner)
        opts = max_inner;
    elseif ~isempty(max_inner)
        opts.max_inner = max_inner;
    end

    opts = fill_method_9_apl_options(opts);

    t_start = tic;
    stats = init_method_9_apl_stats(opts.max_inner);
    beta = opts.beta;
    theta = opts.theta;

    x = helper_project_simplex(x(:), D);
    b = b(:);

    fbest = norm(A * x - b, 1);
    LB = affine_lb_prob1_pd(A, b, x, D);
    xbest = x;

    if fbest - LB < eps
        stats.exit_reason = 'initial_eps';
        stats.cpu_time = toc(t_start);
        fprintf('Initial point reached eps: Ub=%.6e, LB=%.6e, cpu=%.3f\n', fbest, LB, toc(t_start));
        return
    end

    for s = 1:S
        [xbest, fbest, LB, phase_info] = APL_reduction_prob1_pd( ...
            xbest, A, b, D, beta, theta, LB, eps, opts.max_inner, t_start, opts.time_limit);
        stats.phases_run = s;
        stats.last_phase_iterations = phase_info.iterations;
        stats.any_hit_inner_cap = stats.any_hit_inner_cap || phase_info.hit_cap;
        stats.num_phases_hit_inner_cap = stats.num_phases_hit_inner_cap + double(phase_info.hit_cap);

        if phase_info.hit_time_limit
            stats.exit_reason = 'time_limit';
            stats.cpu_time = toc(t_start);
            fprintf('ratio = %.4f\n', (fbest - LB) / max(eps, 1e-16));
            disp('----out of time----');
            return
        end

        if fbest - LB < eps
            stats.exit_reason = 'eps';
            stats.cpu_time = toc(t_start);
            fprintf('Reached eps: s=%d, Ub=%.6e, LB=%.6e, cpu=%.3f\n', s, fbest, LB, toc(t_start));
            return
        end

        if toc(t_start) > opts.time_limit
            stats.exit_reason = 'time_limit';
            stats.cpu_time = toc(t_start);
            fprintf('ratio = %.4f\n', (fbest - LB) / max(eps, 1e-16));
            disp('----out of time----');
            return
        end
    end

    stats.exit_reason = 'completed_S';
    stats.cpu_time = toc(t_start);
end

function [xupper, Ub, LB, info] = APL_reduction_prob1_pd(x, A, b, D, beta, theta, LB_in, eps, max_inner, t_start, time_limit)
    k = 1;

    x = helper_project_simplex(x(:), D);
    center = x;
    xupper = x;
    LB0 = LB_in;

    Ub0 = norm(A * xupper - b, 1);
    Ub = Ub0;
    LB = LB_in;

    l = beta * LB + (1 - beta) * Ub0;
    threshold = (1 - theta) * l + theta * Ub0;
    while Ub > threshold && k <= max_inner
        if toc(t_start) > time_limit
            info.iterations = k - 1;
            info.hit_cap = false;
            info.hit_time_limit = true;
            return
        end

        xupper_prev = xupper;
        Ub_prev = Ub;

        alpha = 2 / (k + 1);
        xlower = (1 - alpha) * xupper + alpha * x;

        r = A * xlower - b;
        g = A' * sign(r);
        fvalue = norm(r, 1) - g' * xlower;
        rhs_level = l - fvalue;

        a_prox = x - center;
        rhs_prox = a_prox' * x;

        [cut_lb, hit_time_limit] = bounded_affine_lb_barx( ...
            g, fvalue, D, a_prox, rhs_prox, t_start, time_limit);
        if hit_time_limit
            info.iterations = k - 1;
            info.hit_cap = false;
            info.hit_time_limit = true;
            return
        end
        LB = max(LB, min(l, cut_lb));

        % Lan's lower-bound exit: stop the phase when the new lower bound
        % reaches l - theta * (l - lb_0), and output the previous upper point.
        if LB >= l - theta * (l - LB0)
            xupper = xupper_prev;
            Ub = Ub_prev;
            info.iterations = k;
            info.hit_cap = false;
            info.hit_time_limit = false;
            return
        end

        [x, is_feasible, hit_time_limit] = proj_level_prox_pd( ...
            center, D, g, rhs_level, a_prox, rhs_prox, t_start, time_limit);
        if hit_time_limit
            info.iterations = k - 1;
            info.hit_cap = false;
            info.hit_time_limit = true;
            return
        end
        if ~is_feasible
            xupper = xupper_prev;
            Ub = Ub_prev;
            info.iterations = k;
            info.hit_cap = false;
            info.hit_time_limit = false;
            return
        end
        x = x(:);

        tol_level = 1e-6 * max(1, abs(rhs_level));
        tol_prox = 1e-6 * max(1, abs(rhs_prox));

        if g' * x > rhs_level + tol_level
            error('method_9_apl:levelCutViolation', ...
                'Projected point violates the level cut.');
        end
        if norm(a_prox) >= 1e-12 && a_prox' * x < rhs_prox - tol_prox
            error('method_9_apl:proxCutViolation', ...
                'Projected point violates the prox cut.');
        end

        xtest = (1 - alpha) * xupper + alpha * x;
        ftest = norm(A * xtest - b, 1);

        if ftest < Ub
            xupper = xtest;
            Ub = ftest;
        end

        if Ub - LB < eps
            fprintf('Inner reached eps: k=%d, Ub=%.6e, LB=%.6e, cpu=%.3f\n', k, Ub, LB, toc(t_start));
            info.iterations = k;
            info.hit_cap = false;
            info.hit_time_limit = false;
            return
        end

        k = k + 1;
    end

    info.iterations = min(k, max_inner);
    info.hit_cap = (k > max_inner) && (Ub > threshold);
    info.hit_time_limit = false;
    if k > max_inner
        warning('method_9_apl:InnerCap', ...
            'method_9_apl inner loop hit max iterations; stopping phase.');
    end
end

function stats = init_method_9_apl_stats(max_inner)
    stats = struct( ...
        'exit_reason', '', ...
        'cpu_time', 0, ...
        'phases_run', 0, ...
        'last_phase_iterations', 0, ...
        'any_hit_inner_cap', false, ...
        'num_phases_hit_inner_cap', 0, ...
        'max_inner', max_inner);
end

function opts = fill_method_9_apl_options(opts)
    defaults.beta = 0.3;
    defaults.theta = 0.2;
    defaults.max_inner = 100;
    defaults.time_limit = 1800;

    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end
end

function lb = affine_lb_prob1_pd(A, b, x, D)
    r = A * x - b;
    g = A' * sign(r);
    fvalue = norm(r, 1) - g' * x;
    lb = bounded_affine_lb(g, fvalue, D);
end

function lb = bounded_affine_lb(g, fvalue, D)
    lb = fvalue + D * min([0; g(:)]);
end

function [lb, hit_time_limit] = bounded_affine_lb_barx(g, fvalue, D, a_prox, rhs_prox, t_start, time_limit)
    hit_time_limit = false;
    if norm(a_prox) < 1e-12
        lb = bounded_affine_lb(g, fvalue, D);
        return
    end

    [sub_lb, hit_time_limit] = linear_lb_simplex_halfspace( ...
        g(:), a_prox(:), rhs_prox, D, t_start, time_limit);
    if hit_time_limit
        lb = NaN;
        return
    end
    lb = fvalue + sub_lb;
end

function [xproj, is_feasible, hit_time_limit] = proj_level_prox_pd( ...
        center, D, g, rhs_level, a_prox, rhs_prox, t_start, time_limit)
    hit_time_limit = false;
    if norm(a_prox) < 1e-12
        [xproj, is_feasible, hit_time_limit] = project_one_cut_simplex( ...
            center, g, rhs_level, D, t_start, time_limit);
        return
    end

    tol_level = 1e-7 * max(1, abs(rhs_level));
    tol_prox = 1e-7 * max(1, abs(rhs_prox));

    x0 = helper_project_simplex(center, D);
    if g' * x0 <= rhs_level + tol_level && a_prox' * x0 >= rhs_prox - tol_prox
        xproj = x0;
        is_feasible = true;
        return
    end

    [x_level, ok_level, hit_time_limit] = project_one_cut_simplex( ...
        center, g, rhs_level, D, t_start, time_limit);
    if hit_time_limit
        xproj = [];
        is_feasible = false;
        return
    end
    if ok_level && a_prox' * x_level >= rhs_prox - tol_prox
        xproj = x_level;
        is_feasible = true;
        return
    end

    [x_prox, ok_prox, hit_time_limit] = project_one_cut_simplex( ...
        center, -a_prox, -rhs_prox, D, t_start, time_limit);
    if hit_time_limit
        xproj = [];
        is_feasible = false;
        return
    end
    if ok_prox && g' * x_prox <= rhs_level + tol_level
        xproj = x_prox;
        is_feasible = true;
        return
    end

    mu_lo = 0;
    if ok_level
        prox_lo = a_prox' * x_level - rhs_prox;
        x_lo = x_level;
    else
        prox_lo = -inf;
        x_lo = [];
    end

    if prox_lo >= -tol_prox
        xproj = x_lo;
        is_feasible = true;
        return
    end

    mu_hi = 1;
    x_hi = [];
    prox_hi = -inf;
    hi_found = false;
    for iter = 1:60
        if toc(t_start) > time_limit
            xproj = [];
            is_feasible = false;
            hit_time_limit = true;
            return
        end
        [x_try, ok_try, hit_time_limit] = project_one_cut_simplex( ...
            center + mu_hi * a_prox, g, rhs_level, D, t_start, time_limit);
        if hit_time_limit
            xproj = [];
            is_feasible = false;
            return
        end
        if ok_try
            prox_hi = a_prox' * x_try - rhs_prox;
            x_hi = x_try;
            if prox_hi >= -tol_prox
                hi_found = true;
                break
            end
        end
        mu_hi = 2 * mu_hi;
    end

    if ~hi_found
        xproj = [];
        is_feasible = false;
        return
    end

    for iter = 1:60
        if toc(t_start) > time_limit
            xproj = [];
            is_feasible = false;
            hit_time_limit = true;
            return
        end
        mu_mid = (mu_lo + mu_hi) / 2;
        [x_mid, ok_mid, hit_time_limit] = project_one_cut_simplex( ...
            center + mu_mid * a_prox, g, rhs_level, D, t_start, time_limit);
        if hit_time_limit
            xproj = [];
            is_feasible = false;
            return
        end
        if ~ok_mid
            mu_lo = mu_mid;
            continue
        end
        prox_mid = a_prox' * x_mid - rhs_prox;
        if prox_mid >= 0
            mu_hi = mu_mid;
            x_hi = x_mid;
            prox_hi = prox_mid;
        else
            mu_lo = mu_mid;
            x_lo = x_mid;
            prox_lo = prox_mid;
        end
    end

    if ~isempty(x_hi) && g' * x_hi <= rhs_level + tol_level && a_prox' * x_hi >= rhs_prox - tol_prox
        xproj = x_hi;
        is_feasible = true;
    else
        xproj = [];
        is_feasible = false;
    end
end

function [xproj, is_feasible, hit_time_limit] = project_one_cut_simplex( ...
        center, c, rhs, D, t_start, time_limit)
    hit_time_limit = false;
    tol = 1e-8 * max(1, abs(rhs));
    x0 = helper_project_simplex(center, D);
    if c' * x0 <= rhs + tol
        xproj = x0;
        is_feasible = true;
        return
    end

    lo = 0;
    hi = 1;
    found = false;
    for iter = 1:80
        if toc(t_start) > time_limit
            xproj = [];
            is_feasible = false;
            hit_time_limit = true;
            return
        end
        xhi = helper_project_simplex(center - hi * c, D);
        if c' * xhi <= rhs + tol
            found = true;
            break
        end
        hi = 2 * hi;
    end

    if ~found
        xproj = [];
        is_feasible = false;
        return
    end

    for iter = 1:80
        if toc(t_start) > time_limit
            xproj = [];
            is_feasible = false;
            hit_time_limit = true;
            return
        end
        mid = (lo + hi) / 2;
        xmid = helper_project_simplex(center - mid * c, D);
        if c' * xmid <= rhs + tol
            hi = mid;
            xproj = xmid;
        else
            lo = mid;
        end
    end

    is_feasible = c' * xproj <= rhs + tol;
end

function [val, hit_time_limit] = linear_lb_simplex_halfspace(g, a, rhs, D, t_start, time_limit)
    hit_time_limit = false;
    tol = 1e-12;
    candidates = 0;

    pos_idx = find(a > tol);
    if ~isempty(pos_idx)
        candidates = [candidates; g(pos_idx) ./ a(pos_idx)];
    end

    n = length(g);
    for i = 1:n
        if toc(t_start) > time_limit
            val = NaN;
            hit_time_limit = true;
            return
        end
        for j = i + 1:n
            denom = a(i) - a(j);
            if abs(denom) > tol
                mu = (g(i) - g(j)) / denom;
                if mu >= 0
                    candidates(end + 1, 1) = mu; %#ok<AGROW>
                end
            end
        end
    end

    candidates = unique(max(candidates, 0));
    if isempty(candidates) || candidates(1) > 0
        candidates = [0; candidates];
    end

    best_dual = -inf;
    for idx = 1:length(candidates)
        if toc(t_start) > time_limit
            val = NaN;
            hit_time_limit = true;
            return
        end
        mu = candidates(idx);
        best_dual = max(best_dual, lb_dual_value(mu, g, a, rhs, D));

        if idx < length(candidates)
            mu_next = candidates(idx + 1);
            mu_mid = 0.5 * (mu + mu_next);
            best_dual = max(best_dual, lb_dual_value(mu_mid, g, a, rhs, D));
        else
            tail_mu = max(mu + 1, 2 * max(1, mu));
            best_dual = max(best_dual, lb_dual_value(tail_mu, g, a, rhs, D));
        end
    end

    if ~isfinite(best_dual)
        error('method_9_apl:lbSubproblemInfeasible', ...
            'Could not solve the lower-bound subproblem without a solver.');
    end

    val = best_dual;
end

function dual_val = lb_dual_value(mu, g, a, rhs, D)
    coeff = g - mu * a;
    dual_val = mu * rhs + D * min([0; coeff(:)]);
end

function tau = level_stepsize(center, g, rhs_level, D)
    tau = 0;
    x0 = helper_project_simplex(center, D);
    if g' * x0 <= rhs_level + 1e-8
        return
    end

    lo = 0;
    hi = 1;
    while g' * helper_project_simplex(center - hi * g, D) > rhs_level + 1e-8
        hi = 2 * hi;
        if hi > 1e12
            break
        end
    end

    if g' * helper_project_simplex(center - hi * g, D) > rhs_level + 1e-8
        error('method_9_apl:infeasibleLevelCut', ...
            'Could not find a feasible point for the level-cut subproblem.');
    end

    for iter = 1:60
        mid = (lo + hi) / 2;
        xmid = helper_project_simplex(center - mid * g, D);
        if g' * xmid <= rhs_level + 1e-8
            hi = mid;
            tau = mid;
        else
            lo = mid;
        end
    end

    xfinal = helper_project_simplex(center - tau * g, D);
    if g' * xfinal > rhs_level + 1e-8
        error('method_9_apl:levelSearchFailed', ...
            'Level-cut line search did not produce a feasible point.');
    end
end

function alpha = pd_two_cut_alpha(center, D, g, rhs_level, a_prox, rhs_prox)
    % Find a convex combination of the level-cut and prox-cut normals whose
    % projected point satisfies both inequalities. This is still a search
    % heuristic, but it now guarantees that the returned alpha corresponds
    % to a feasible projected point.
    xg = helper_project_simplex(center - g, D);
    xa = helper_project_simplex(center - a_prox, D);

    fg = a_prox' * xg - rhs_prox;
    fa = g' * xa - rhs_level;
    lg = g' * xg - rhs_level;
    pa = a_prox' * xa - rhs_prox;

    if lg <= 1e-8 && fg >= -1e-8
        alpha = 1;
        return
    end

    if fa <= 1e-8 && pa >= -1e-8
        alpha = 0;
        return
    end

    left = 0;
    right = 1;
    alpha = 0.5;
    for iter = 1:40
        alpha = (left + right) / 2;
        xmid = helper_project_simplex(center - (alpha * g + (1 - alpha) * a_prox), D);
        gap = g' * xmid - rhs_level;

        if gap > 0
            left = alpha;
        else
            right = alpha;
        end
    end

    best_alpha = NaN;
    best_dist = inf;
    alpha_grid = linspace(0, 1, 81);
    candidates = unique([alpha_grid, alpha, right, left, 1, 0]);

    for alpha_try = candidates
        xtry = helper_project_simplex(center - (alpha_try * g + (1 - alpha_try) * a_prox), D);
        if g' * xtry <= rhs_level + 1e-8 && a_prox' * xtry >= rhs_prox - 1e-8
            dist = norm(xtry - center);
            if dist < best_dist
                best_dist = dist;
                best_alpha = alpha_try;
            end
        end
    end

    if ~isnan(best_alpha)
        alpha = best_alpha;
        return
    end

    error('method_9_apl:twoCutSearchFailed', ...
        'Could not find a feasible projected point for the level/prox subproblem.');
end

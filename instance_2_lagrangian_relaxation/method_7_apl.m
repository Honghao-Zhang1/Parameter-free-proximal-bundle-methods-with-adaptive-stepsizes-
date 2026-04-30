function [pi, fbest, LB, stats] = method_7_apl(q, W, T, h, x, pi_0, opt, S, opts)
    % Accelerated proximal level method for the Lagrangian dual in prob2.
    %
    % Sign convention:
    % - helper_dual_oracle returns the original dual objective d(pi), which we want
    %   to maximize and compare directly with the known primal optimum opt.
    % - Internally, APL is run on Ub(pi) = -d(pi), a convex minimization
    %   problem with exact lower bound -opt.
    % - On output, fbest is converted back to the original dual value and
    %   LB is returned as the original known primal optimum opt.
    %
    % To make the expensive MIP oracle practical while staying close to APL,
    % the implementation uses:
    % - exact lower-bound information,
    % - more aggressive default level parameters,
    % - early exits when an inner phase stops improving.

    if nargin < 9 || isempty(opts)
        opts = struct();
    end
    opts = fill_method_7_apl_options(opts);

    tic;
    stats = init_method_7_apl_stats(opts.max_inner);

    x = x(:);
    pi = pi_0(:);
    opt_primal = opt;
    LB_internal = -opt_primal;

    [dual_value, g] = helper_dual_oracle(q, W, T, h, pi, x);
    fbest_internal = -dual_value;

    if fbest_internal - LB_internal < opts.eps
        stats.exit_reason = 'initial_optimal';
        stats.cpu_time = toc;
        fbest = dual_value;
        LB = opt_primal;
        fprintf('Initial point optimal: dual=%.6e\n', fbest);
        return
    end

    stalled_phases = 0;
    hit_inner_cap = false;
    for s = 1:S
        gap_before = fbest_internal - LB_internal;

        [pi, fbest_internal, g, phase_info] = APL_reduction_lar( ...
            q, W, T, h, x, pi, fbest_internal, g, LB_internal, opts);
        hit_inner_cap = hit_inner_cap || phase_info.hit_cap;
        stats.phases_run = s;
        stats.any_hit_inner_cap = hit_inner_cap;
        stats.num_phases_hit_inner_cap = stats.num_phases_hit_inner_cap + double(phase_info.hit_cap);
        stats.last_phase_iterations = phase_info.iterations;
        stats.last_phase_rel_improve = phase_info.rel_improve;

        if fbest_internal - LB_internal < opts.eps
            stats.exit_reason = 'eps';
            stats.cpu_time = toc;
            fbest = -fbest_internal;
            LB = opt_primal;
            if hit_inner_cap
                fprintf('method_7_apl note: at least one phase hit max_inner=%d.\n', opts.max_inner);
            end
            fprintf('Reached optimality: s=%d, dual=%.6e, opt=%.6e\n', s, fbest, LB);
            return
        end

        if phase_info.rel_improve <= opts.phase_stall_tol
            stalled_phases = stalled_phases + 1;
            if stalled_phases >= opts.phase_stall_limit
                stats.exit_reason = 'phase_stall';
                stats.cpu_time = toc;
                fprintf(['method_7_apl stopped early after %d stalled phases: ' ...
                    'gap=%.6e, last phase iters=%d, cpu=%.3f\n'], ...
                    stalled_phases, opt_primal + fbest_internal, phase_info.iterations, toc);
                fbest = -fbest_internal;
                LB = opt_primal;
                if hit_inner_cap
                    fprintf('method_7_apl note: at least one phase hit max_inner=%d.\n', opts.max_inner);
                end
                return
            end
        else
            stalled_phases = 0;
        end

        if toc > opts.time_limit
            stats.exit_reason = 'time_limit';
            stats.cpu_time = toc;
            fprintf('ratio = %.4f\n', (fbest_internal - LB_internal) / max(opts.eps, 1e-16));
            disp('----out of time----');
            fbest = -fbest_internal;
            LB = opt_primal;
            if hit_inner_cap
                fprintf('method_7_apl note: at least one phase hit max_inner=%d.\n', opts.max_inner);
            end
            return
        end

        % Guard against pathological phases that appear to improve only in
        % floating-point noise.
        if gap_before - (fbest_internal - LB_internal) <= opts.abs_phase_progress
            stalled_phases = stalled_phases + 1;
            if stalled_phases >= opts.phase_stall_limit
                stats.exit_reason = 'absolute_progress_stall';
                stats.cpu_time = toc;
                fprintf('method_7_apl stopped: phase progress below absolute threshold.\n');
                fbest = -fbest_internal;
                LB = opt_primal;
                if hit_inner_cap
                    fprintf('method_7_apl note: at least one phase hit max_inner=%d.\n', opts.max_inner);
                end
                return
            end
        end
    end

    stats.exit_reason = 'completed_S';
    stats.cpu_time = toc;
    fbest = -fbest_internal;
    LB = opt_primal;
    if hit_inner_cap
        fprintf('method_7_apl note: at least one phase hit max_inner=%d.\n', opts.max_inner);
    end
end

function stats = init_method_7_apl_stats(max_inner)
    stats = struct( ...
        'exit_reason', '', ...
        'cpu_time', 0, ...
        'phases_run', 0, ...
        'last_phase_iterations', 0, ...
        'last_phase_rel_improve', 0, ...
        'any_hit_inner_cap', false, ...
        'num_phases_hit_inner_cap', 0, ...
        'max_inner', max_inner);
end

function [pi_upper, Ub, g_upper, info] = APL_reduction_lar(q, W, T, h, x, pi, Ub, g_upper, LB, opts)
    k = 1;
    center = pi;
    pi_upper = pi;
    Ub0 = Ub;

    level = opts.beta * LB + (1 - opts.beta) * Ub0;
    threshold = (1 - opts.theta) * level + opts.theta * Ub0;

    no_improve = 0;
    while Ub > threshold && k <= opts.max_inner
        alpha = 2 / (k + 1);
        pi_lower = (1 - alpha) * pi_upper + alpha * pi;

        [fmax_lower, g_lower] = helper_dual_oracle(q, W, T, h, pi_lower, x);
        r_lower = -fmax_lower;

        rhs_level = level - (r_lower + g_lower' * pi_lower);
        a_prox = pi - center;
        rhs_prox = a_prox' * pi;

        if norm(a_prox) < opts.prox_tol
            pi_candidate = proj_onecut_lar(center, -g_lower, rhs_level, 'le', opts.proj_tol);
        else
            pi_candidate = proj_twocuts_lar(center, -g_lower, rhs_level, a_prox, rhs_prox, opts.proj_tol);
        end

        pi_test = (1 - alpha) * pi_upper + alpha * pi_candidate;
        [fmax_test, g_test] = helper_dual_oracle(q, W, T, h, pi_test, x);
        ftest = -fmax_test;

        if ftest < Ub - opts.abs_phase_progress
            pi_upper = pi_test;
            Ub = ftest;
            g_upper = g_test;
            no_improve = 0;
        else
            no_improve = no_improve + 1;
        end

        if Ub - LB < opts.eps
            fprintf('Inner reached eps: k=%d, Ub=%.6e, cpu=%.3f\n', k, Ub, toc);
            break
        end

        pi = pi_candidate;

        if no_improve >= opts.inner_stall_limit
            break
        end

        k = k + 1;
    end

    info.iterations = k;
    info.hit_cap = (k > opts.max_inner) && (Ub > threshold);
    info.rel_improve = max(0, (Ub0 - Ub) / max(1, abs(Ub0 - LB)));
end

function opts = fill_method_7_apl_options(opts)
    defaults.eps = 1e-4;
    defaults.beta = 0.9;
    defaults.theta = 0.35;
    defaults.max_inner = 40;
    defaults.time_limit = 1800;
    defaults.inner_stall_limit = 5;
    defaults.phase_stall_limit = 3;
    defaults.phase_stall_tol = 1e-3;
    defaults.abs_phase_progress = 1e-8;
    defaults.proj_tol = 1e-8;
    defaults.prox_tol = 1e-12;
    defaults.norm_tol = 1e-14;

    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end
end

function pi = proj_onecut_lar(center, c, d, sense, tol)
    % Project onto a single affine halfspace:
    %   If sense='le': c'*z <= d
    %   If sense='ge': c'*z >= d

    pi = center;
    val = c' * pi;
    c_norm_sq = c' * c;

    if c_norm_sq <= tol^2
        return
    end

    if strcmpi(sense, 'le')
        if val <= d + tol
            return
        end
        pi = center - ((val - d) / c_norm_sq) * c;
    else
        if val >= d - tol
            return
        end
        pi = center + ((d - val) / c_norm_sq) * c;
    end
end

function pi = proj_twocuts_lar(center, g, r_level, a, b_prox, tol)
    % Project onto the two halfspaces
    %   g'*z <= r_level
    %   a'*z >= b_prox

    pi = center;
    if g' * pi <= r_level + tol && a' * pi >= b_prox - tol
        return
    end

    pi_g = proj_onecut_lar(center, g, r_level, 'le', tol);
    if a' * pi_g >= b_prox - tol
        pi = pi_g;
        return
    end

    pi_a = proj_onecut_lar(center, a, b_prox, 'ge', tol);
    if g' * pi_a <= r_level + tol
        pi = pi_a;
        return
    end

    gg = g' * g;
    aa = a' * a;
    ga = g' * a;
    detM = gg * aa - ga * ga;

    % For an unconstrained domain, the exact projection onto the
    % intersection of two halfspaces is obtained by an active-set solve:
    % - center itself if feasible,
    % - projection onto one active halfspace if that already satisfies the
    %   other,
    % - otherwise the projection onto both active boundaries.
    if detM >= 1e-18 * max(1, gg * aa)
        rhs1 = g' * center - r_level;
        rhs2 = b_prox - a' * center;

        lam = (aa * rhs1 + ga * rhs2) / detM;
        mu = (ga * rhs1 + gg * rhs2) / detM;
        pi = center - lam * g + mu * a;

        if lam >= -tol && mu >= -tol && g' * pi <= r_level + tol && a' * pi >= b_prox - tol
            return
        end
    end

    % Degenerate case: the two normals are (nearly) collinear, so the
    % feasible set reduces to a slab or a single halfspace along g.
    if gg <= tol^2
        pi = pi_a;
        return
    end

    kappa = ga / gg;
    if abs(kappa) <= tol
        pi = pi_g;
        return
    end

    if kappa > 0
        lower = b_prox / kappa;
        upper = r_level;
        t = g' * center;

        if lower <= upper + tol
            if t < lower - tol
                pi = center + ((lower - t) / gg) * g;
            elseif t > upper + tol
                pi = center - ((t - upper) / gg) * g;
            else
                pi = center;
            end
        elseif norm(pi_g - center) <= norm(pi_a - center)
            pi = pi_g;
        else
            pi = pi_a;
        end
    else
        upper = min(r_level, b_prox / kappa);
        pi = proj_onecut_lar(center, g, upper, 'le', tol);
    end
end

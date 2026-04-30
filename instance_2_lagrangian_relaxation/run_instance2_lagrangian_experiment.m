function [results, data] = run_instance2_lagrangian_experiment(opts)
%RUN_INSTANCE2_LAGRANGIAN_EXPERIMENT Reproduce the Lagrangian-relaxation comparison experiment.
%
%   [results, data] = run_instance2_lagrangian_experiment() generates the
%   random instance from the legacy Instance 2 script, builds the master
%   iterates, and runs the requested comparison methods on subproblems from
%   iterates 2:4.
%
%   Optional opts fields:
%
%       opts.seed             Random seed. Default: 15.
%       opts.n                Number of x/y/u variables. Default: 120.
%       opts.m_master         Number of master constraints. Default: 50.
%       opts.m_lar            Number of Lagrangian constraints. Default: 5.
%       opts.num_subproblems  Number of q vectors. Default: 20.
%       opts.master_iters     Number of master iterations. Default: 4.
%       opts.eval_iters       Master iterates to test. Default: 2:4.
%       opts.methods          Methods to run. Default: all comparison methods.
%       opts.tau              Adaptive bundle contraction parameter. Default: 0.95.
%       opts.eps              Target accuracy for subgradient/PB. Default: 1e-4.
%       opts.beta             PB sufficient decrease parameter. Default: 0.1.
%       opts.lam_min          PB minimum stepsize. Default: 1e-6.
%       opts.a                PB growth/shrink factor. Default: 2.
%       opts.K_pb             PB iteration cap. Default: 10000.
%       opts.S_apl            APL phase cap. Default: 10000.
%       opts.apl_opts         Options passed to method_7_apl.
%
%   Method names are:
%       'subgradient', 'Ad-GPB', 'P-Ad-GPB', 'PB', 'APL'.

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = fill_instance2_lagrangian_options(opts);

    data = generate_instance2_lagrangian_data(opts);
    fprintf('Instance 2: n=%d, q=%d, master_iters=%d\n', ...
        opts.n, opts.num_subproblems, opts.master_iters);

    [data.x_iterate, data.opt_val, data.Bundle_pi, data.Bundle_L] = ...
        build_master_iterates(data, opts);

    results = struct();
    results.options = opts;
    results.master = struct( ...
        'x_iterate', data.x_iterate, ...
        'opt_val', data.opt_val, ...
        'Bundle_pi', data.Bundle_pi, ...
        'Bundle_L', data.Bundle_L);

    if should_run(opts, 'subgradient')
        disp('---------subgradient-----------')
        results.subgradient = run_subgradient_comparison(data, opts);
        disp(['cpu=', num2str(results.subgradient.cpu_time)])
    end

    if should_run(opts, 'adgpb')
        disp('---------Ad-GPB-----------')
        results.adgpb = run_pairwise_method(data, opts, ...
            @(q, x, opt) method_3_adgpb(q, data.W, data.T, data.h, x, opt, data.pi_0, opts.tau));
        disp(['cpu=', num2str(results.adgpb.cpu_time)])
    end

    if should_run(opts, 'p_adgpb')
        disp('---------P-Ad-GPB-----------')
        results.p_adgpb = run_pairwise_method(data, opts, ...
            @(q, x, opt) method_4_p_adgpb(q, data.W, data.T, data.h, x, opt, data.pi_0, opts.tau, []));
        disp(['cpu=', num2str(results.p_adgpb.cpu_time)])
    end

    if should_run(opts, 'pb')
        disp('---------PB-----------')
        results.pb = run_pairwise_method(data, opts, ...
            @(q, x, opt) method_5_pb(data.pi_0, opts.beta, opts.eps, opts.lam_min, ...
                opts.a, q, data.W, data.T, data.h, x, opt, opts.K_pb));
        disp(['cpu=', num2str(results.pb.cpu_time)])
    end

    if should_run(opts, 'apl')
        disp('---------APL-----------')
        results.apl = run_apl_comparison(data, opts);
        disp(['cpu=', num2str(results.apl.cpu_time)])
    end
end

function data = generate_instance2_lagrangian_data(opts)
    rng(opts.seed)

    data.A = randi(100, opts.m_master, opts.n);
    data.C = randi(100, opts.m_master, opts.n);
    data.c = randi(100, opts.n, 1);
    data.T = randi(100, opts.m_lar, opts.n);
    data.W = randi(100, opts.m_lar, opts.n);
    data.d = randi(100, opts.n, 1);
    data.b = 3/4 * (data.A * ones(opts.n, 1) + data.C * ones(opts.n, 1));
    data.h = 3/4 * (data.T * ones(opts.n, 1) + data.W * ones(opts.n, 1));
    data.q = cell(opts.num_subproblems, 1);

    for i = 1:opts.num_subproblems
        data.q{i} = randi(100, opts.n, 1);
    end

    data.x_0 = ones(opts.n, 1);
    data.pi_0 = rand(opts.n, 1);
end

function [x_iterate, opt_val, Bundle_pi, Bundle_L] = build_master_iterates(data, opts)
    x = data.x_0;
    x_iterate = zeros(opts.master_iters, opts.n);
    opt_val = zeros(opts.num_subproblems * opts.master_iters, 1);
    Bundle_pi = zeros(opts.master_iters, opts.n);
    Bundle_L = zeros(opts.master_iters, 1);

    t_start = tic;
    for k = 1:opts.master_iters
        x_iterate(k, :) = x;
        for i = 1:opts.num_subproblems
            subopt = helper_primal_reference(data.q{i}, data.W, data.T, data.h, x);
            opt_val(subproblem_index(k, i, opts), 1) = subopt;
            [pi, cut_value] = method_4_p_adgpb( ...
                data.q{i}, data.W, data.T, data.h, x, subopt, data.pi_0, opts.tau, []);
            Bundle_pi(k, :) = Bundle_pi(k, :) + pi';
            Bundle_L(k, 1) = Bundle_L(k, 1) + cut_value;
        end

        Bundle_pi(k, :) = Bundle_pi(k, :) / opts.num_subproblems;
        Bundle_L(k, 1) = Bundle_L(k, 1) / opts.num_subproblems;
        x = helper_original_problem(data.c, data.d, data.A, data.C, data.b, Bundle_pi, Bundle_L, k);
    end
    fprintf('Master generation cpu=%g\n', toc(t_start));
end

function stats = run_subgradient_comparison(data, opts)
    t_start = tic;
    calls = 0;

    for k = opts.eval_iters
        x = data.x_iterate(k, :)';
        for j = 1:opts.num_subproblems
            pi = data.pi_0;
            opt = data.opt_val(subproblem_index(k, j, opts));
            [fpi, g] = helper_dual_oracle(data.q{j}, data.W, data.T, data.h, pi, x);
            calls = calls + 1;

            while opt - fpi > opts.eps
                pi = pi + (opt - fpi) / norm(g)^2 * g;
                [fpi, g] = helper_dual_oracle(data.q{j}, data.W, data.T, data.h, pi, x);
                calls = calls + 1;
            end
        end
    end

    stats = struct('cpu_time', toc(t_start), 'oracle_calls', calls);
end

function stats = run_pairwise_method(data, opts, method_handle)
    t_start = tic;
    calls = 0;

    for k = opts.eval_iters
        x = data.x_iterate(k, :)';
        for j = 1:opts.num_subproblems
            opt = data.opt_val(subproblem_index(k, j, opts));
            method_handle(data.q{j}, x, opt);
            calls = calls + 1;
        end
    end

    stats = struct('cpu_time', toc(t_start), 'subproblem_runs', calls);
end

function stats = run_apl_comparison(data, opts)
    t_start = tic;
    calls = 0;
    last_stats = struct();

    for k = opts.eval_iters
        x = data.x_iterate(k, :)';
        for j = 1:opts.num_subproblems
            opt = data.opt_val(subproblem_index(k, j, opts));
            [~, ~, ~, last_stats] = method_7_apl( ...
                data.q{j}, data.W, data.T, data.h, x, data.pi_0, opt, opts.S_apl, opts.apl_opts);
            calls = calls + 1;
        end
    end

    stats = struct( ...
        'cpu_time', toc(t_start), ...
        'subproblem_runs', calls, ...
        'last_method_stats', last_stats);
end

function idx = subproblem_index(k, j, opts)
    idx = opts.num_subproblems * (k - 1) + j;
end

function tf = should_run(opts, method)
    requested = cellfun(@canonical_method_name, opts.methods, 'UniformOutput', false);
    tf = any(strcmpi(requested, canonical_method_name(method)));
end

function name = canonical_method_name(name)
    name = char(name);
    key = lower(strrep(strrep(strrep(name, '-', ''), '*', ''), '_', ''));

    switch key
        case {'sub', 'subgradient', 'polyaksubgradient'}
            name = 'subgradient';
        case {'adgpb', 'method3', 'method3adgpb'}
            name = 'adgpb';
        case {'padgpb', 'poladgpb', 'method4', 'method4padgpb'}
            name = 'p_adgpb';
        case {'pb', 'gcpb', 'gcpblar', 'method5', 'method5pb'}
            name = 'pb';
        case {'apl', 'method7', 'method7apl'}
            name = 'apl';
        otherwise
            name = key;
    end
end

function opts = fill_instance2_lagrangian_options(opts)
    defaults.seed = 15;
    defaults.n = 120;
    defaults.m_master = 50;
    defaults.m_lar = 5;
    defaults.num_subproblems = 20;
    defaults.master_iters = 4;
    defaults.eval_iters = 2:4;
    defaults.methods = {'subgradient', 'Ad-GPB', 'P-Ad-GPB', 'PB', 'APL'};
    defaults.tau = 0.95;
    defaults.eps = 1e-4;
    defaults.beta = 0.1;
    defaults.lam_min = 1e-6;
    defaults.a = 2;
    defaults.K_pb = 10000;
    defaults.S_apl = 10000;
    defaults.apl_opts = struct();

    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end

    if opts.n ~= 120 || opts.m_master ~= 50
        error('run_instance2_lagrangian_experiment:unsupportedSize', ...
            'helper_original_problem currently supports only n=120 and m_master=50.');
    end

    if ischar(opts.methods) || isstring(opts.methods)
        opts.methods = cellstr(opts.methods);
    end
end

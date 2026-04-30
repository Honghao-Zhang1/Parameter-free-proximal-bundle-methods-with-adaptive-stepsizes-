function [results, data] = run_instance1_sparse_l1_experiment(d, n, s, opts)
%RUN_INSTANCE1_SPARSE_L1_EXPERIMENT Generate and solve sparse l1 instances.
%
%   [results, data] = run_instance1_sparse_l1_experiment(d, n, s) creates
%   a reproducible sparse least-absolute-deviation instance
%
%       min ||A*x - b||_1  subject to x >= 0,
%
%   then runs the comparison methods used for instance 1.
%
%   Optional opts fields:
%
%       opts.seed          Random seed. Default: 11.
%       opts.time_limit    Time limit for Pol-Sub. Default: 3600.
%       opts.K_bundle      Iteration cap for GPB/Ad-GPB methods. Default: 100000.
%       opts.K_classical   Iteration cap for PB/BG methods. Default: 10000.
%       opts.S_apl         Number of APL outer phases. Default: 10000.
%       opts.beta          Classical bundle sufficient decrease parameter. Default: 0.1.
%       opts.tau           Adaptive GPB contraction parameter. Default: 0.95.
%       opts.lam_min       Minimum stepsize for PB. Default: 1e-6.
%       opts.a             PB growth/shrink factor. Default: 2.
%       opts.pol_c         Polynomial adaptive GPB scale. Default: 40.
%       opts.apl_opts      Options passed to method_7_apl.
%       opts.methods       Cell array of method names to run. Default: all.
%
%   Method names are:
%       'Pol-Sub', 'Ad-GPB*', 'GPB', 'P-Ad-GPB*',
%       'PB', 'GCPBBG', 'APL'.
%
%   This driver builds the sparse nonnegative l1 feasibility instance
%   A = D*N, where N is sparse Gaussian and D is a random diagonal row
%   scaling. The stopping target is relative to the initial objective value.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = fill_sparse_l1_experiment_options(opts);

    data = generate_sparse_prob1(d, n, s, opts.seed);
    A = data.A;
    b = data.b;
    x_0 = data.x0;
    x_true = data.x_true;

    fprintf('d=%d, n=%d, s=%g\n', d, n, s);

    eps_target = 1e-4 * norm(A*x_0' - b, 1);
    initial_lam = norm(A*x_0' - b, 1) / norm(A' * sign(A*x_0' - b))^2;

    results = struct();
    results.eps = eps_target;
    results.initial_lam = initial_lam;
    results.methods = opts.methods;

    if should_run(opts, 'subgradient')
        disp('----------Pol-Sub----------')
        results.subgradient = method_1_polyak_subgradient(x_0, eps_target, A, b, opts.time_limit);
    end

    if should_run(opts, 'adgpb')
        disp('----------Ad-GPB*----------')
        results.adgpb = run_printed_method(@() method_3_adgpb( ...
            opts.K_bundle, initial_lam, x_0, eps_target, opts.tau, n, A, b));
    end

    if should_run(opts, 'gpb')
        disp('----------GPB----------')
        results.gpb = run_printed_method(@() method_2_gpb( ...
            opts.K_bundle, initial_lam, x_0, eps_target, n, A, b));
    end

    if should_run(opts, 'pol_adgpb')
        disp('----------P-Ad-GPB*----------')
        results.pol_adgpb1 = run_printed_method(@() method_4_p_adgpb_star( ...
            opts.K_bundle, x_0, eps_target, opts.tau, n, A, b, opts.pol_c));
    end

    if should_run(opts, 'pb')
        disp('----------PB*----------')
        results.pb = run_printed_method(@() method_5_pb( ...
            initial_lam, x_0, opts.beta, eps_target, A, b, ...
            opts.K_classical, n, opts.lam_min, opts.a));
    end

    if should_run(opts, 'gcpbbg')
        disp('----------GCPBBG----------')
        results.gcpbbg = run_printed_method(@() method_6_bg( ...
            x_0, opts.beta, eps_target, A, b, opts.K_classical, n, x_true));
    end

    if should_run(opts, 'apl')
        disp('----------APL----------')
        t_start = tic;
        [x_apl, apl_stats] = method_7_apl(x_0, eps_target, A, b, opts.S_apl, opts.apl_opts);
        results.apl = struct( ...
            'x', x_apl, ...
            'fvalue', norm(A*x_apl - b, 1), ...
            'stats', apl_stats, ...
            'cpu_time', toc(t_start));
    end
end

function data = generate_sparse_prob1(d, n, s, seed)
%GENERATE_SPARSE_PROB1 Build the reproducible random test instance.

    rng(seed)
    A = sprandn(d, n, s);

    % Row scaling makes the nonsmooth objective less uniformly conditioned
    % than the raw sparse Gaussian matrix.
    row_scale = 1000 * diag(rand(d, 1));
    A = sparse(row_scale * A);

    x_true = randn(1, n).^2;
    b = A * x_true';
    x0 = rand(1, n).^2;

    data = struct('A', A, 'b', b, 'x0', x0, 'x_true', x_true);
end

function stats = run_printed_method(method_handle)
%RUN_PRINTED_METHOD Time legacy methods that report progress with disp.

    t_start = tic;
    method_handle();
    stats = struct('cpu_time', toc(t_start));
end

function tf = should_run(opts, method)
    requested = cellfun(@canonical_method_name, opts.methods, 'UniformOutput', false);
    tf = any(strcmpi(requested, canonical_method_name(method)));
end

function name = canonical_method_name(name)
    name = char(name);
    key = lower(strrep(strrep(strrep(name, '-', ''), '*', ''), '_', ''));

    switch key
        case {'polsub', 'subgradient', 'polyaksubgradient'}
            name = 'subgradient';
        case 'adgpb'
            name = 'adgpb';
        case 'gpb'
            name = 'gpb';
        case {'padgpb', 'poladgpb', 'poladgpb1'}
            name = 'pol_adgpb';
        case {'pb', 'gcpb', 'gcpb2'}
            name = 'pb';
        case 'gcpbbg'
            name = 'gcpbbg';
        case 'apl'
            name = 'apl';
        otherwise
            name = key;
    end
end

function opts = fill_sparse_l1_experiment_options(opts)
    defaults.seed = 11;
    defaults.time_limit = 3600;
    defaults.K_bundle = 100000;
    defaults.K_classical = 10000;
    defaults.S_apl = 10000;
    defaults.beta = 0.1;
    defaults.tau = 0.95;
    defaults.lam_min = 1e-6;
    defaults.a = 2;
    defaults.pol_c = 40;
    defaults.apl_opts = struct();
    defaults.methods = { ...
        'Pol-Sub', 'Ad-GPB*', 'GPB', 'P-Ad-GPB*', ...
        'PB', 'GCPBBG', 'APL'};

    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end

    if ischar(opts.methods) || isstring(opts.methods)
        opts.methods = cellstr(opts.methods);
    end
end

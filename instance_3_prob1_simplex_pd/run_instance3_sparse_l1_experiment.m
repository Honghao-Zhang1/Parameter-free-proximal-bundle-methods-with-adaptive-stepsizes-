function [xbest, fbest, LB, stats] = run_instance3_sparse_l1_experiment(d,n,s,D,c,opts)
%RUN_INSTANCE3_SPARSE_L1_EXPERIMENT Generate and solve a sparse prob1 benchmark instance.
%
%   [xbest, fbest, LB, stats] = run_instance3_sparse_l1_experiment(d, n, s, D, c)
%   creates a reproducible sparse matrix A using sprandn(d,n,s), scales the
%   rows, and runs comparison methods on
%
%       min ||A*x - b||_1  subject to x >= 0, sum(x) <= D.
%
%   c is kept as an input for compatibility with older experiment calls. The
%   active methods in this file do not use c.
%
%   Optional opts fields:
%
%       opts.seed        Random seed. Default: 11.
%       opts.K_bundle    Iteration cap for Ad-GPB. Default: 30000.
%       opts.S           Number of APL outer phases. Default: 10000.
%       opts.tau         Ad-GPB contraction parameter. Default: 0.95.
%       opts.apl_opts    Options passed directly to method_9_apl.
%       opts.methods     Cell array of method names. Default: {'Ad-GPB','APL'}.

if nargin < 6 || isempty(opts)
    opts = struct();
end
opts = fill_instance3_options(opts);

rng(opts.seed)
disp(['d=',num2str(d), ', n = ', num2str(n), ', s = ', num2str(s),', D=', num2str(D), ', c=', num2str(c)])

% Generate a sparse, row-scaled linear model. The row scaling creates a more
% ill-conditioned nonsmooth objective than a plain sprandn instance.
A = sprandn(d,n,s);
row_scale = 1000 * diag(rand(d,1));
A = sparse(row_scale * A);

b = randn(d,1).^2;
x_0 = 1/D*ones(1,n);




%% Problem scale and target accuracy
g =  A'*sign(A*x_0'-b);
fx = norm(A*x_0'-b,1);

% Linearize ||A*x-b||_1 at x_0 and minimize that affine model over
% {x >= 0, sum(x) <= D}; this gives the initial lower bound used to scale
% the requested accuracy.
if all(g>=0)
    Lb = fx-g'*x_0' ;
else
    Lb = fx-g'*x_0' + D*min(g);
end
eps = 1e-3* (fx-Lb) ;
initial_lam = (fx-Lb)/(norm(g)^2);

stats = struct();
stats.eps = eps;
stats.initial_lam = initial_lam;
stats.initial_ub = fx;
stats.initial_lb = Lb;
stats.methods = opts.methods;

%% Ad-GPB
if should_run(opts, 'adgpb')
    disp('----------Ad-GPB----------')
    t_start = tic;
    method_8_adgpb(opts.K_bundle, initial_lam, x_0, eps, opts.tau, n, A, b, D);
    stats.adgpb = struct('cpu_time', toc(t_start));
end

%% APL method
if should_run(opts, 'apl')
    disp('----------APL----------')
    [xbest, fbest, LB, apl_stats] = method_9_apl(x_0, eps, A, b, opts.S, D, opts.apl_opts);
    stats.apl = apl_stats;
    fprintf('Final gap = %.6e, exit = %s\n', fbest - LB, apl_stats.exit_reason);
else
    xbest = [];
    fbest = [];
    LB = [];
end
end

function tf = should_run(opts, method)
    requested = cellfun(@canonical_method_name, opts.methods, 'UniformOutput', false);
    tf = any(strcmpi(requested, canonical_method_name(method)));
end

function name = canonical_method_name(name)
    name = char(name);
    key = lower(strrep(strrep(strrep(name, '-', ''), '*', ''), '_', ''));

    switch key
        case {'adgpb', 'method8', 'method8adgpb'}
            name = 'adgpb';
        case {'apl', 'method9', 'method9apl'}
            name = 'apl';
        otherwise
            name = key;
    end
end

function opts = fill_instance3_options(opts)
    defaults.seed = 11;
    defaults.K_bundle = 30000;
    defaults.S = 10000;
    defaults.tau = 0.95;
    defaults.apl_opts = struct();
    defaults.methods = {'Ad-GPB', 'APL'};

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

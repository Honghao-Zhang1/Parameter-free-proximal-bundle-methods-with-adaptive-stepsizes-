function x = helper_project_nonnegative_two_cuts(x0,g,r,a,b,tol)
%HELPER_PROJECT_NONNEGATIVE_TWO_CUTS Project onto R_+^n with two affine cuts.
%
%   x = helper_project_nonnegative_two_cuts(x0, g, r, a, b, tol) computes the Euclidean
%   projection of x0 onto the intersection of the nonnegative orthant, the
%   level cut g'*x <= r, and the prox/aggregate cut a'*x >= b.
%
%   This helper supports the two-cut bundle projection specialized to
%   h = indicator{x >= 0}. It is a subproblem solver, not a standalone
%   method.

if nargin<6, tol=1e-10; end
x0=x0(:); g=g(:); a=a(:);
proj=@(z) max(z,0);

% try no cuts
x = proj(x0);
if is_feasible_point(x, g, r, a, b, tol), return; end

% try only g-cut (1D bisection)
xg = proj_onecut(x0,g,r,'le',tol);
if is_feasible_point(xg, g, r, a, b, tol), x=xg; return; end

% try only a-cut (1D bisection)
xa = proj_onecut(x0,a,b,'ge',tol);
if is_feasible_point(xa, g, r, a, b, tol), x=xa; return; end

% Try the exact both-active KKT system on the current positive support.
[x, ok] = proj_two_active_support(x0, g, r, a, b, tol);
if ok && is_feasible_point(x, g, r, a, b, tol)
    return;
end

% Robust fallback: Dykstra-style projection onto
%   {x>=0} \cap {g'*x <= r} \cap {a'*x >= b}
[x, ok] = proj_two_cuts_dykstra(x0, g, r, a, b, tol);
if ok && is_feasible_point(x, g, r, a, b, tol)
    return;
end

% Exact fallback for rare hard cases.
[x, ok] = proj_two_cuts_mosek(x0, g, r, a, b);
if ok
    return;
end

if norm(xg - x0) <= norm(xa - x0)
    x = xg;
else
    x = xa;
end
end

function x = proj_onecut(x0,c,d,sense,tol)
x0=x0(:); c=c(:); proj=@(z) max(z,0);
x=proj(x0); val=c'*x;
if strcmpi(sense,'le')
    if val<=d+tol, return; end
    lo=0; hi=1;
    while c'*proj(x0-hi*c) > d+tol, hi=2*hi; end
    for k=1:60
        mid=(lo+hi)/2;
        xm=proj(x0-mid*c);
        if c'*xm <= d+tol, hi=mid; x=xm; else, lo=mid; end
    end
else
    if val>=d-tol, return; end
    lo=0; hi=1;
    while c'*proj(x0+hi*c) < d-tol, hi=2*hi; end
    for k=1:20
        mid=(lo+hi)/2;
        xm=proj(x0+mid*c);
        if c'*xm >= d-tol, hi=mid; x=xm; else, lo=mid; end
    end
end
end

function x = proj_halfspace(x0,c,d,sense,tol)
c = c(:);
cn2 = c' * c;
if cn2 <= tol^2
    x = x0;
    return;
end
val = c' * x0;
if strcmpi(sense,'le')
    if val <= d + tol
        x = x0;
    else
        x = x0 - ((val - d) / cn2) * c;
    end
else
    if val >= d - tol
        x = x0;
    else
        x = x0 + ((d - val) / cn2) * c;
    end
end
end

function [x, ok] = proj_two_active_support(x0, g, r, a, b, tol)
x = [];
ok = false;
active = x0 > tol;

for it = 1:50
    idx = find(active);
    if isempty(idx)
        return;
    end

    gS = g(idx);
    aS = a(idx);
    xS = x0(idx);

    gg = gS' * gS;
    aa = aS' * aS;
    ga = gS' * aS;
    gx = gS' * xS;
    ax = aS' * xS;

    M = [gg, -ga; ga, -aa];
    rhs = [gx - r; ax - b];
    if rcond(M) <= 1e-12
        return;
    end

    sol = M \ rhs;
    lam = sol(1);
    mu = sol(2);
    if lam < -1e-9 || mu < -1e-9
        return;
    end

    x = max(x0 - lam * g + mu * a, 0);
    new_active = x > 10 * tol;

    if isequal(new_active, active)
        eq1 = abs(g' * x - r);
        eq2 = abs(a' * x - b);
        if is_feasible_point(x, g, r, a, b, tol) && ...
                eq1 <= 1e-7 * max(1, abs(r)) && ...
                eq2 <= 1e-7 * max(1, abs(b))
            ok = true;
        end
        return;
    end

    active = new_active;
end
end

function [x, ok] = proj_two_cuts_dykstra(x0, g, r, a, b, tol)
x = x0;
ok = false;
p1 = zeros(size(x0));
p2 = zeros(size(x0));
p3 = zeros(size(x0));
best_feas = [];
best_dist = Inf;

for it=1:2000
    x_old = x;

    y = max(x + p1, 0);
    p1 = x + p1 - y;

    z = proj_halfspace(y + p2, g, r, 'le', tol);
    p2 = y + p2 - z;

    x = proj_halfspace(z + p3, a, b, 'ge', tol);
    p3 = z + p3 - x;

    xcand = max(x, 0);
    if is_feasible_point(xcand, g, r, a, b, tol)
        dist = norm(xcand - x0);
        if dist < best_dist
            best_feas = xcand;
            best_dist = dist;
        end
    end

    if norm(x - x_old) <= 1e-10 * max(1, norm(x_old)) && ~isempty(best_feas)
        x = best_feas;
        ok = true;
        return;
    end
end

if ~isempty(best_feas)
    x = best_feas;
    ok = true;
end
end

function [x, ok] = proj_two_cuts_mosek(x0, g, r, a, b)
ok = false;
x = [];

if exist('mosekopt', 'file') ~= 2
    mosek_path = 'C:/Program Files/Mosek/10.0/toolbox/r2017a';
    if exist(mosek_path, 'dir') == 7
        addpath(mosek_path);
    end
end
if exist('mosekopt', 'file') ~= 2
    return;
end

n = numel(x0);
clear prob;
prob.blx = zeros(n, 1);
prob.bux = [];
prob.c = -x0(:);
prob.qosubi = (1:n)';
prob.qosubj = (1:n)';
prob.qoval = ones(n, 1);
prob.a = [g(:)'; -a(:)'];
prob.blc = [-inf; -inf];
prob.buc = [r; -b];

[rcode, res] = mosekopt('minimize echo(0)', prob);
if rcode ~= 0 || ~isfield(res, 'sol') || ~isfield(res.sol, 'itr') || ...
        ~isfield(res.sol.itr, 'xx') || isempty(res.sol.itr.xx)
    return
end

if isfield(res.sol.itr, 'solsta')
    solsta = res.sol.itr.solsta;
    if strcmp(solsta, 'PRIMAL_INFEASIBLE') || strcmp(solsta, 'PRIMAL_INFEASIBLE_CER')
        return
    end
    if ~(strcmp(solsta, 'OPTIMAL') || strcmp(solsta, 'NEAR_OPTIMAL'))
        return
    end
end

x = res.sol.itr.xx(:);
ok = true;
end

function tf = is_feasible_point(x, g, r, a, b, tol)
tf = all(x >= -10 * tol) && ...
    g' * x <= r + 1e-7 * max(1, abs(r)) + tol && ...
    a' * x >= b - 1e-7 * max(1, abs(b)) - tol;
end

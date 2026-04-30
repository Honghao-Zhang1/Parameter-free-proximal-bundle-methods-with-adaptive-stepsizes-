function xP = helper_project_simplex(x,D)
%helper_project_simplex Project a vector onto {x >= 0, sum(x) <= D}.
%
%   xP = helper_project_simplex(x, D) returns the Euclidean projection of x onto the
%   nonnegative l1 ball with radius D. If x is already feasible, it is
%   returned unchanged. Otherwise the projection is computed by sorting and
%   finding the simplex threshold.

    input_size = size(x);
    x_col = x(:);
    n = length(x_col);

    if sum(x_col) <= D && all(x_col >= 0)
        xP = x;
        return
    end

    u = sort(x_col, 'descend');
    csU = cumsum(u);
    iPos = (u + (D - csU) ./ (1:n)') > 0;
    rho = find(iPos, 1, 'last');
    if isempty(rho)
        xP = zeros(input_size);
        return
    end

    lambda = (D - csU(rho)) / rho;
    xP = reshape(max(x_col + lambda, 0), input_size);
end

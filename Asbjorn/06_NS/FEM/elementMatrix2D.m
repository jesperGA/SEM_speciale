function [A, B, grad1, grad2, J] = elementMatrix2D(xy, xi, w, D, N)
    % elementMatrix2D: Calculate 2D element matrices A and B
    %
    % Arguments:
    % x, y (arrays): Coordinates.
    % xi (array): Grid points.
    % w (array): Weights.
    % D (matrix): Derivative matrix.
    % N (int): Polynomial order.
    %
    % Returns:
    % A, B (matrices): 2D element matrices.

    x = reshape(xy(:, 1), N + 1, N + 1);
    y = reshape(xy(:, 2), N + 1, N + 1);

    ldof_u = (length(x)) ^ 2;

    A_4D = zeros(N + 1, N + 1, N + 1, N + 1);
    B_4D = zeros(N + 1, N + 1, N + 1, N + 1);
    A = zeros(ldof_u, ldof_u);
    B = zeros(ldof_u, ldof_u);

    grad1 = zeros(ldof_u, ldof_u);
    grad2 = zeros(ldof_u, ldof_u);

    [xr, xs, yr, ys] = calculateDerivatives(x, y, D, N);
    J = xr .* ys - xs .* yr;

    grad = calculateGradient(xr, xs, yr, ys, D, N);

    for i = 1:N + 1
        for j = 1:N + 1
            for m = 1:N + 1
                for n = 1:N + 1
                    A_4D(i, j, m, n) = computeAElement(i, j, m, n, N, w, J, grad);
                    B_4D(i, j, m, n) = computeBElement(i, j, m, n, N, w, J);

                    row = (i - 1) * (N + 1) + j;
                    col = (m - 1) * (N + 1) + n;

                    grad1(row, col) = 1/J(i,j)*grad(i, j, m, n, 1);
                    grad2(row, col) = 1/J(i,j)*grad(i, j, m, n, 2);

                    A(row, col) = A_4D(i, j, m, n);
                    B(row, col) = B_4D(i, j, m, n);
                end
            end
        end
    end
end

function [xr, xs, yr, ys] = calculateDerivatives(x, y, D, N)
    % Calculate derivatives of x and y coordinates
    xr = zeros(N + 1, N + 1);
    xs = zeros(N + 1, N + 1);
    yr = zeros(N + 1, N + 1);
    ys = zeros(N + 1, N + 1);

    for p = 1:N + 1
        for q = 1:N + 1
            for m = 1:N + 1
                xr(p, q) = xr(p, q) + D(p, m) * x(m, q);
                xs(p, q) = xs(p, q) + D(q, m) * x(p, m);
                yr(p, q) = yr(p, q) + D(p, m) * y(m, q);
                ys(p, q) = ys(p, q) + D(q, m) * y(p, m);
            end
        end
    end
end

function grad = calculateGradient(xr, xs, yr, ys, D, N)
    % Calculate gradient
    grad = zeros(N + 1, N + 1, N + 1, N + 1, 2);
    for p = 1:N + 1
        for q = 1:N + 1
            for m = 1:N + 1
                for n = 1:N + 1
                    d_qn = q == n; % Kronecker delta
                    d_pm = p == m; % Kronecker delta
                    grad(p,q,m,n,1) = xr(p,q) * d_pm * D(q,n) - xs(p,q) * D(p,m) * d_qn;
                    grad(p,q,m,n,2) = ys(p,q) * D(p,m) * d_qn - yr(p,q) * d_pm * D(q,n);
                end
            end
        end
    end
end

function aElem = computeAElement(i, j, m, n, N, w, J, grad)
    % Compute A matrix element
    aElem = 0;
    for p = 1:N + 1
        for q = 1:N + 1
            aElem = aElem + w(p) * w(q) * 1 / J(p, q) * ...
                    dot(grad(p, q, i, j, :), grad(p, q, m, n, :));
        end
    end
end

function bElem = computeBElement(i, j, m, n, N, w, J)
    % Compute B matrix element
    d_im = i == m; % Kronecker delta
    d_jn = j == n; % Kronecker delta
    bElem = w(i) * w(j) * abs(J(i, j)) * d_im * d_jn;
end

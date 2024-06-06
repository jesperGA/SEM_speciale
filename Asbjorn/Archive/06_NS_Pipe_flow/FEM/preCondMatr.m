function [M] = preCondMatr(w, N, xy, ldof);
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

    L1 = max(xy(:,1)) - min(xy(:,1));
    L2 = max(xy(:,2)) - min(xy(:,2));

    M_4D = zeros(N + 1, N + 1, N + 1, N + 1);
    M = zeros(ldof, ldof);

    J = (L1/2)*(L2/2);

    for i = 1:N + 1
        for j = 1:N + 1
            for m = 1:N + 1
                for n = 1:N + 1
                    M_4D(i, j, m, n) = computeBElement(i, j, m, n, N, w, J);

                    row = (i - 1) * (N + 1) + j;
                    col = (m - 1) * (N + 1) + n;
                    M(row, col) = M_4D(i, j, m, n);
                end
            end
        end
    end
end

function bElem = computeBElement(i, j, m, n, N, w, J)
    % Compute B matrix element
    d_im = i == m; % Kronecker delta
    d_jn = j == n; % Kronecker delta
    bElem = w(i) * w(j) * abs(J) * d_im * d_jn;
end
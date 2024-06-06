function [fit] = lagrange2D(X, Y, N, nodalValues, XI, YI)
    % N = 5; % Order of the polynomial (degree N-1)
    % 
    % [X, Y] = meshgrid(GetGLL(N+1), GetGLL(N+1)); % LGL nodes in 2D
    % 
    % % Define a function to interpolate (for demonstration)
    % F = @(x, y) sin(pi*x).*cos(pi*y); % Example function
    % nodalValues = F(X, Y); % Evaluate function at nodal points
    % nodalValues =randi([0 1], N+1,N+1); % Evaluate function at nodal points
    % 
    % % Points where we want to interpolate
    % [XI, YI] = meshgrid(linspace(0, 1, 6), linspace(0, 1, 6));
    % 
    % Perform the interpolation
    fit = zeros(size(XI));
    for ix = 1:size(XI, 1)
        for iy = 1:size(XI, 2)
            fit(ix, iy) = interpolateValue(N, XI(ix, iy), YI(ix, iy), X, Y, nodalValues);
        end
    end
end

function value = interpolateValue(N, xi, yi, X, Y, nodalValues)
    % Interpolate the value at (xi, yi) using Lagrange polynomials
    Lx = lagrangePolynomials(xi, X(1,:));
    Ly = lagrangePolynomials(yi, Y(:,1));
    
    value = 0;
    for i = 1:N+1
        for j = 1:N+1
            value = value + nodalValues(i,j) * Lx(i) * Ly(j);
        end
    end
end

function L = lagrangePolynomials(x, nodes)
    % Compute Lagrange polynomials for all nodes at point x
    n = length(nodes);
    L = zeros(1, n);
    for i = 1:n
        L(i) = 1;
        for j = 1:n
            if i ~= j
                L(i) = L(i) * (x - nodes(j)) / (nodes(i) - nodes(j));
            end
        end
    end
end

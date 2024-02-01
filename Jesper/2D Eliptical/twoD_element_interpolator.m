function [interp_data] = twoD_element_interpolator(mesh,mesh_interp, solution, X, Y)
% Assuming solution is now a matrix with rows corresponding to elements
% and columns to solution values at each node of the element
IX = mesh.IX;
interp_data = [];
IX_interp = mesh_interp.IX;
xG = mesh_interp.X(:,2);
yG = mesh_interp.X(:,3);
for e = 1:size(IX,3)
    % Extract node indices for element 'e'
    nen = IX(:,:,e);
    n = length(nen);

    % Get coordinates for these nodes
    x_data = X(nen);
    y_data = Y(nen);
    x_data = x_data(:);
    y_data = y_data(:);



    nen_grid = IX_interp(:,:,e);
    xg = xG(nen_grid);
    yg = yG(nen_grid);
    
    n_grid = length(nen_grid);

    x_grid = reshape(xg,n_grid,n_grid);
    y_grid = reshape(yg,n_grid,n_grid);


    % Define 2D grid for interpolation

    % [x_grid, y_grid] = meshgrid(linspace(min(x_data), max(x_data)), ...
    %                             linspace(min(y_data'), max(y_data')));

    % Perform 2D interpolation
    z_data = reshape(solution(nen, :),n,n); % Solution values at element nodes
    z_fit = Lagrange2DInterpolation(x_data, y_data', z_data, x_grid, y_grid);

    % Store interpolated values
    interp_data{e} = [x_grid(:), y_grid(:), z_fit(:)];
end
end

function z_fit = Lagrange2DInterpolation(x_data, y_data, z_data, x_grid, y_grid)
% Check if the number of data points in x, y and z are consistent
% if length(x_data) ~= length(y_data) || length(y_data) ~= size(z_data, 1)
%     disp('The number of data points must be consistent in x, y, and z.');
%     z_fit = NaN;
%     return;
% end

% Number of data points
n = length(x_data);

% Initialize output matrix
z_fit = zeros(size(x_grid));

% Loop over each grid point
for gx = 1:size(x_grid, 1)
    for gy = 1:size(x_grid, 2)
        % Current grid point
        x = x_grid(gx, gy);
        y = y_grid(gx, gy);

        % Compute the 2D Lagrange polynomial for the current grid point
        z_value = 0;
        for i = 1:n
            for j = 1:n
                % Compute the Lagrange basis polynomial for x_data(i) and y_data(j)
                Lx = 1;
                Ly = 1;
                for k = [1:i-1, i+1:n]
                    Lx = Lx * (x - x_data(k)) / (x_data(i) - x_data(k));
                end
                for k = [1:j-1, j+1:n]
                    Ly = Ly * (y - y_data(k)) / (y_data(j) - y_data(k));
                end

                % Multiply the basis polynomials by the corresponding z value
                z_value = z_value + z_data(i, j) * Lx * Ly;
            end
        end

        % Assign the computed value to the output matrix
        z_fit(gx, gy) = z_value;
    end
end
end

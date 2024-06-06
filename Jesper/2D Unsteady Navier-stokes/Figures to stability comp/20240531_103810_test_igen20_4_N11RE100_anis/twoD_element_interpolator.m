function [interp_data] = twoD_element_interpolator(mesh_data, N, u,v,p)
% [interp_data] = twoD_element_interpolator(IX_data,X_data, mesh_interp, solution)
% Writes datafile that can be intepreted by writenek with interpolated data.
% Assuming solution is now a matrix with rows corresponding to elements
% and columns to solution values at each node of the element
IXv = mesh_data.IXv;
X_data = mesh_data.Xv;
Xp_data = mesh_data.Xp;
IXp = mesh_data.IXp;


interp_data = zeros([size(IXv,3),N^2,5]);
temp_results = cell(size(IXv,3), 1);

parfor e = 1:size(IXv,3)
    % Extract node indices for element 'e'

    nenV = IXv(:,:,e);
    nV = length(nenV);

    % Get coordinates for these nodes
    x_data = X_data(nenV,2);
    y_data = X_data(nenV,3);
    x_data = reshape(x_data,nV,nV);
    y_data = reshape(y_data,nV,nV);

    nenP = IXp(:,:,e);
    nP = length(nenP);

    % Get coordinates for these nodes
    xp_data = Xp_data(nenP,2);
    yp_data = Xp_data(nenP,3);
    xp_data = reshape(xp_data,nP,nP);
    yp_data = reshape(yp_data,nP,nP);

    xmax = max(X_data(nenV,2));
    xmin = min(X_data(nenV,2));

    ymax = max(X_data(nenV,3));
    ymin = min(X_data(nenV,3));

    xg = linspace(xmin,xmax,N);
    yg = linspace(ymin,ymax,N);

    [x_grid,y_grid] = meshgrid(xg,yg);

    % Perform 2D interpolation
    u_data = reshape(u(nenV, :),nV,nV); % Solution values at element nodes
    uu = Lagrange2DInterpolation(x_data, y_data', u_data, x_grid, y_grid);

    v_data = reshape(v(nenV, :),nV,nV); % Solution values at element nodes
    vv = Lagrange2DInterpolation(x_data, y_data', v_data, x_grid, y_grid);

    % Perform 2D interpolation
    p_data = reshape(p(nenP, :),nP,nP); % Solution values at element nodes
    pp = Lagrange2DInterpolation(xp_data, yp_data', p_data, x_grid, y_grid);

    % Store interpolated values
    % interp_data(e,:,1) = x_grid(:)';
    % interp_data(e,:,2) = y_grid(:)';
    % interp_data(e,:,3) = uu(:)';
    % interp_data(e,:,4) = vv(:)';
    % interp_data(e,:,5) = pp(:)';
    % interp_data{e} = [x_grid(:), y_grid(:), z_fit(:)];

    temp_results{e} = [x_grid(:), y_grid(:), uu(:), vv(:), pp(:)];
end
for e = 1:size(IXv, 3)
    % Ensure the size of the matrix from the cell matches the target size
    current_data = temp_results{e};  % Fetch the matrix from cell

    % Assign each slice of the third dimension individually
    interp_data(e, :, 1) = current_data(:, 1)';  % x_grid values
    interp_data(e, :, 2) = current_data(:, 2)';  % y_grid values
    interp_data(e, :, 3) = current_data(:, 3)';  % uu values
    interp_data(e, :, 4) = current_data(:, 4)';  % vv values
    interp_data(e, :, 5) = current_data(:, 5)';  % pp values


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

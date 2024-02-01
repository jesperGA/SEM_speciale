function [interp_data] = twoD_element_interpolator(mesh, solution, X, Y)
    % Assuming solution is now a matrix with rows corresponding to elements
    % and columns to solution values at each node of the element
    IX = mesh.IX; 
    interp_data = [];
    for e = 1:size(IX,3)
        % Extract node indices for element 'e'
        nen = IX(:,:,e);
        
        % Get coordinates for these nodes
        x_data = X(nen, 2);
        y_data = Y(nen, 2);

        % Define 2D grid for interpolation
        [x_grid, y_grid] = meshgrid(linspace(min(x_data), max(x_data), n_grid_pts), ...
                                    linspace(min(y_data), max(y_data), n_grid_pts));

        % Perform 2D interpolation
        z_data = solution(e, :); % Solution values at element nodes
        z_fit = Lagrange2DInterpolation(x_data, y_data, z_data, x_grid, y_grid);

        % Store interpolated values
        interp_data{e} = [x_grid(:), y_grid(:), z_fit(:)];
    end
end
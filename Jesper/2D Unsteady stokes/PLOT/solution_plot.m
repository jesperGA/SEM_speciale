function solution_plot(interp_data,n)

% Assume 'interp_data' is the output from twoD_element_interpolator
figure();
hold on

for e = 1:length(interp_data)
    % Extract data for the e-th element
    data = interp_data{e};
    x_grid = reshape(data(:, 1), [n, n]);
    y_grid = reshape(data(:, 2), [n, n]);
    z_fit = reshape(data(:, 3), [n, n]);

    % Plot using surf
    
    surf(x_grid, y_grid, z_fit);
    
end

xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Interpolated Value');
    title(['Interpolation of Element ', num2str(e)]);
    colorbar;
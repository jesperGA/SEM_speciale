% Parameters
N = 15; % Change this to set grid resolution
R = 0.4; % Radius from the center point

% Generate grid points
x = linspace(0, 1, N);
y = linspace(0, 1, N);
[X, Y] = meshgrid(x, y);

% Calculate the distance from the center point (0.5, 0.5)
distances = sqrt((X - 0.5).^2 + (Y - 0.5).^2);

% Create a mask where distances are less than R
mask = distances < R;

% Plotting
figure;
hold on;
axis equal;
% title('Colored Grid with Circle at Radius R from Center');

% Iterate over each cell in the grid
for i = 1:N-1
    for j = 1:N-1
        % Define the corners of the current cell
        x_corners = [x(i) x(i+1) x(i+1) x(i)];
        y_corners = [y(j) y(j) y(j+1) y(j+1)];
        
        % Check if the current cell is within the radius
        if mask(i,j) || mask(i+1,j) || mask(i,j+1) || mask(i+1,j+1)
            % Color the cell if any of its corners are within the radius
            fill(x_corners, y_corners, 'r');
        else
            % Optionally draw a line around cells outside the radius
            plot(x_corners([1:4 1]), y_corners([1:4 1]), 'k');
        end
    end
end

% Draw the actual circle outline
theta = linspace(0, 2*pi, 200);
x_circle = 0.5 + R * cos(theta);
y_circle = 0.5 + R * sin(theta);
plot(x_circle, y_circle, 'k', 'LineWidth', 2); % Blue circle outline

% Remove the axis for a cleaner look
axis off;
hold off;

function plotmesh(iglob,x,y)

% Example MATLAB script to plot x and y coordinates

% Assuming iglob, x, and y are already defined in your workspace

% Determine the size of iglob
[dim1, dim2, dim3] = size(iglob);

% Initialize arrays to store x and y coordinates
all_x = [];
all_y = [];

% Loop through each element in iglob to extract x and y coordinates
for e = 1:dim3
    for i = 1:dim1
        for j = 1:dim2
            % Get the node index
            nodeIndex = iglob(i, j, e);
            
            % Extract x and y coordinates
            node_x = x(nodeIndex);
            node_y = y(nodeIndex);
            
            % Append to the arrays
            all_x = [all_x, node_x];
            all_y = [all_y, node_y];
        end
    end
end

% Plotting the coordinates
scatter(all_x, all_y);
title('Scatter Plot of x and y Coordinates');
xlabel('x Coordinate');
ylabel('y Coordinate');
grid on;
end
function [N] = find_normal(nodes,X)
% [N] = find_normal(nodes, X)
%
% This function computes the normal vectors for the specified boundary nodes.
% The function assumes that the lines are defined by the boundary nodes and
% calculates normals that are consistent in direction based on the node ordering.
%
% Assistance in generating basic logic in this code was provided by ChatGPT, a large
% language model developed by OpenAI. 
%
%
% Input:
%   nodes - A vector of indices for the nodes on the boundary of the mesh.
%           These nodes should be along the set of lines for which normals are needed.
%   X     - An N x 3 matrix where each row represents a node;
%           the first column contains node indices,
%           the second and third columns contain the x and y coordinates of the nodes, respectively.
%
% Output:
%   N     - An M x 2*M matrix containing the normal vectors for each boundary node specified in 'nodes',
%           where M is the length of 'nodes'. Each row in N corresponds to the normal vector
%           of the boundary node at the same position in 'nodes'. Matrix is
%           used for the calculation N*V, that computes the normal
%           direction of the velocity field.
neqn = size(X,1);
N = zeros(neqn,2*neqn);
n1 = zeros(neqn,1);n2 = n1;

nodes = unique(nodes);

domain_nodes = setdiff(X(:,1),nodes);
X_inside = X(domain_nodes,:);
for i = 1:length(nodes)
    % Extract the coordinates of node i
    node_i_coords = X(nodes(i), 2:3);

    % Compute the Euclidean distances from node i to all other nodes
    distances = sqrt(sum((X(nodes, 2:3) - node_i_coords).^2, 2));

    % Set the distance from node i to itself to infinity to avoid selecting it as the closest node
    distances(distances==0) = inf;

    % Find the index of the closest node
    [~, sub_index] = min(distances);
    closest_node_index = nodes(sub_index);

    % Coordinates of the first node (e.g., node i)
    x1 = X(nodes(i), 2);
    y1 = X(nodes(i), 3);

    % Coordinates of the second node (e.g., the closest node)
    x2 = X(closest_node_index, 2);
    y2 = X(closest_node_index, 3);

    % Compute the direction vector from node i to the closest node
    direction_vector = [x2 - x1, y2 - y1];

    % Compute the normal vector by rotating the direction vector by 90 degrees
    normal_vector = [direction_vector(2), -direction_vector(1)];
    normal_vector = normal_vector/norm(normal_vector,2);
    %%
    %Find correct sign by defining which way into the domain.
    distances = sqrt(sum((X_inside(:, 2:3) - node_i_coords).^2, 2));
    distances(distances==0) = inf;
    % Find the index of the closest node
    [~, sub_index] = min(distances);
    closest_node_index = X_inside(sub_index,1);
    % Coordinates of the first node (e.g., node i)
    x1 = X(nodes(i), 2);
    y1 = X(nodes(i), 3);

    % Coordinates of the second node (e.g., the closest node)
    x2 = X(closest_node_index, 2);
    y2 = X(closest_node_index, 3);

    % Compute the direction vector from node i to the closest node
    direction_vector = [x2 - x1, y2 - y1];
    N_sign = sign(direction_vector);
    
    %Check if the vectors point in the same direction
    if dot(normal_vector,N_sign) > 0
    
        normal_vector = -normal_vector;

    end

    %for corner nodes, check colinearty. 
    cp = normal_vector(1)*direction_vector(2)-normal_vector(2)*direction_vector(1);
    if abs(cp)>1e-6
        normal_vector = -direction_vector/norm(direction_vector,2);
    end


    % normal_vector = N_sign.*normal_vector;

    N(nodes(i), [nodes(i),neqn+nodes(i)]) = normal_vector;

    n1(nodes(i)) = normal_vector(1); n2(nodes(i)) = normal_vector(2);

end

end
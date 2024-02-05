function mesh_plot(mesh)

figure;
hold on;

for e = 1:size(mesh.IXv, 3) % For each element in 'v' mesh
    % Plot interior grid lines with thinner lines
    % Loop through all but the outer rows and columns for the internal grid
    for i = 2:size(mesh.IXv, 1)-1 % For each interior column
        for j = 1:size(mesh.IXv, 2)-1 % Connect nodes vertically
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i, j+1, e);
            plot([mesh.Xv(nodeStart, 2), mesh.Xv(nodeEnd, 2)], [mesh.Xv(nodeStart, 3), mesh.Xv(nodeEnd, 3)], 'k--', 'LineWidth', 0.5);
        end
    end
    for j = 2:size(mesh.IXv, 2)-1 % For each interior row
        for i = 1:size(mesh.IXv, 1)-1 % Connect nodes horizontally
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i+1, j, e);
            plot([mesh.Xv(nodeStart, 2), mesh.Xv(nodeEnd, 2)], [mesh.Xv(nodeStart, 3), mesh.Xv(nodeEnd, 3)], 'k--', 'LineWidth', 0.5);
        end
    end

    % Plot edge lines with thicker lines
    % Top and bottom edge rows
    for i = [1, size(mesh.IXv, 1)] % Only first and last row for edges
        for j = 1:size(mesh.IXv, 2)-1 % Connect nodes horizontally
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i, j+1, e);
            plot([mesh.Xv(nodeStart, 2), mesh.Xv(nodeEnd, 2)], [mesh.Xv(nodeStart, 3), mesh.Xv(nodeEnd, 3)], 'k-', 'LineWidth', 2); % Thicker line
        end
    end
    % Left and right edge columns
    for j = [1, size(mesh.IXv, 2)] % Only first and last column for edges
        for i = 1:size(mesh.IXv, 1)-1 % Connect nodes vertically
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i+1, j, e);
            plot([mesh.Xv(nodeStart, 2), mesh.Xv(nodeEnd, 2)], [mesh.Xv(nodeStart, 3), mesh.Xv(nodeEnd, 3)], 'k-', 'LineWidth', 2); % Thicker line
        end
    end
end

% Plot 'p' mesh with dotted lines
for e = 1:size(mesh.IXp, 3) % For each element in 'p' mesh
    % Assuming a similar logic applies to 'p' mesh, if not, adjust accordingly
    for i = 1:size(mesh.IXp, 1) % For each column
        for j = 1:size(mesh.IXp, 2) - 1 % Connect nodes vertically
            nodeStart = mesh.IXp(i, j, e);
            nodeEnd = mesh.IXp(i, j+1, e);
            plot([mesh.Xp(nodeStart, 2), mesh.Xp(nodeEnd, 2)], [mesh.Xp(nodeStart, 3), mesh.Xp(nodeEnd, 3)], 'k:', 'LineWidth', 0.5);
        end
    end
    for j = 1:size(mesh.IXp, 2) % For each row
        for i = 1:size(mesh.IXp, 1) - 1 % Connect nodes horizontally
            nodeStart = mesh.IXp(i, j, e);
            nodeEnd = mesh.IXp(i+1, j, e);
            plot([mesh.Xp(nodeStart, 2), mesh.Xp(nodeEnd, 2)], [mesh.Xp(nodeStart, 3), mesh.Xp(nodeEnd, 3)], 'k:', 'LineWidth', 0.5);
        end
    end
end
scatter(mesh.Xv(:,2),mesh.Xv(:,3),'*k')
scatter(mesh.Xp(:,2),mesh.Xp(:,3),'or')
hold off;
axis equal;
xlabel('X coordinate');
ylabel('Y coordinate');
title('Staggered Mesh Visualization');
end
function mesh_plot(mesh)
figure;
hold on;

% Initialize arrays for 'p' mesh lines (only interior)
pLinesX = [];
pLinesY = [];

% Initialize arrays to hold line coordinates for efficiency
edgeLinesX = [];
edgeLinesY = [];
interiorLinesX = [];
interiorLinesY = [];

% Process 'v' mesh
for e = 1:size(mesh.IXv, 3) % For each element in 'v' mesh
    % Edge lines
    edgeIdx = [1, size(mesh.IXv,1); 1, size(mesh.IXv,2)]; % Indices for edges (first and last row/column)
    for i = edgeIdx(1,:) % Horizontal edge lines (top and bottom)
        for j = 1:size(mesh.IXv,2)-1
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i, j+1, e);
            edgeLinesX = [edgeLinesX, [mesh.Xv(nodeStart, 2); mesh.Xv(nodeEnd, 2)]];
            edgeLinesY = [edgeLinesY, [mesh.Xv(nodeStart, 3); mesh.Xv(nodeEnd, 3)]];
        end
    end
    for j = edgeIdx(2,:) % Vertical edge lines (left and right)
        for i = 1:size(mesh.IXv,1)-1
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i+1, j, e);
            edgeLinesX = [edgeLinesX, [mesh.Xv(nodeStart, 2); mesh.Xv(nodeEnd, 2)]];
            edgeLinesY = [edgeLinesY, [mesh.Xv(nodeStart, 3); mesh.Xv(nodeEnd, 3)]];
        end
    end

    % Interior lines
    for i = 2:size(mesh.IXv,1)-1 % Horizontal interior lines
        for j = 1:size(mesh.IXv,2)-1
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i, j+1, e);
            interiorLinesX = [interiorLinesX, [mesh.Xv(nodeStart, 2); mesh.Xv(nodeEnd, 2)]];
            interiorLinesY = [interiorLinesY, [mesh.Xv(nodeStart, 3); mesh.Xv(nodeEnd, 3)]];
        end
    end
    for j = 2:size(mesh.IXv,2)-1 % Vertical interior lines
        for i = 1:size(mesh.IXv,1)-1
            nodeStart = mesh.IXv(i, j, e);
            nodeEnd = mesh.IXv(i+1, j, e);
            interiorLinesX = [interiorLinesX, [mesh.Xv(nodeStart, 2); mesh.Xv(nodeEnd, 2)]];
            interiorLinesY = [interiorLinesY, [mesh.Xv(nodeStart, 3); mesh.Xv(nodeEnd, 3)]];
        end
    end
end

% Plot lines for 'v' mesh - edges and interior
plot(edgeLinesX, edgeLinesY, 'k-', 'LineWidth', 2,'HandleVisibility', 'off'); % Edge lines
plot(interiorLinesX, interiorLinesY, 'k-', 'LineWidth', 0.5,'HandleVisibility', 'off'); % Interior lines

% Process 'p' mesh for interior lines
for e = 1:size(mesh.IXp, 3) % For each element in 'p' mesh
    for i = 1:size(mesh.IXp,1) % Horizontal lines
        for j = 1:size(mesh.IXp,2)-1
            nodeStart = mesh.IXp(i, j, e);
            nodeEnd = mesh.IXp(i, j+1, e);
            pLinesX = [pLinesX, [mesh.Xp(nodeStart, 2); mesh.Xp(nodeEnd, 2)]];
            pLinesY = [pLinesY, [mesh.Xp(nodeStart, 3); mesh.Xp(nodeEnd, 3)]];
        end
    end
    for j = 1:size(mesh.IXp,2) % Vertical lines
        for i = 1:size(mesh.IXp,1)-1
            nodeStart = mesh.IXp(i, j, e);
            nodeEnd = mesh.IXp(i+1, j, e);
            pLinesX = [pLinesX, [mesh.Xp(nodeStart, 2); mesh.Xp(nodeEnd, 2)]];
            pLinesY = [pLinesY, [mesh.Xp(nodeStart, 3); mesh.Xp(nodeEnd, 3)]];
        end
    end
end

% Plot lines for 'p' mesh - interior only, using dotted lines
plot(pLinesX, pLinesY, 'k:', 'LineWidth',1,'HandleVisibility', 'off');

scatter(mesh.Xv(:,2),mesh.Xv(:,3),'*k','DisplayName','$v$-grid')
scatter(mesh.Xp(:,2),mesh.Xp(:,3),'or','DisplayName','$p$-grid')
ax = gca; % Get the handle to the current axes
ax.TickLabelInterpreter = 'latex'; % Set the tick labels to use LaTeX interpreter
ax.FontSize = 14; % Set the font size for the axes
legend('Location','northoutside','Orientation','horizontal','Interpreter','latex',FontSize=18)
hold off;
axis equal;
xlabel('X coordinate','Interpreter','latex',FontSize=18);
ylabel('Y coordinate','Interpreter','latex',FontSize=18);
% title('Staggered Mesh Visualization');
end
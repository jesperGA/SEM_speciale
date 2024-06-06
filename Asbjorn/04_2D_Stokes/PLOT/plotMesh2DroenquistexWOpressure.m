clear;
close all;
clc;
warning('off','all');

% Add necessary paths
addpath('FEM');
addpath('MESH'); 
addpath('PLOT');
addpath('MISC');

LX = 2;
LY = LX;
NELX = 1; % Number of elements along each side
NELY = NELX;

study.example = 'Roenquist';

for N = [1 2 3 5]
    study.N = N;
    [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);
    plot_quadrature(mesh, study.N);
end
% 
% study.N = 2;
% plot_quadrature(mesh, study.N);


function plot_quadrature(mesh, N)

fntsize = 24;
linethk = 3;
mrkSize = 14;
defaultColors = get(groot, 'DefaultAxesColorOrder');

n_gl=N+1;
[GL_points, ~]  = lgwt(n_gl,-1,1);

fig = figure;
fig.Position = [651.4000 419.4000 432 430.4000];
hold on;

% Labels and title
xlabel('\it x_1');
ylabel('\it x_2');
axis equal;
grid on;
constant = 0.3;
xlim([min(mesh.X(:,2)) - constant, max(mesh.X(:,2)) + constant]);
ylim([min(mesh.X(:,3)) - constant, max(mesh.X(:,3)) + constant]);

% Initialize arrays to hold line coordinates for efficiency
edgeLinesX = [];
edgeLinesY = [];

% Process 'v' mesh
for e = 1:size(mesh.IX, 3) % For each element in 'v' mesh

    % Edge lines
    edgeIdx = [1, size(mesh.IX, 1); 1, size(mesh.IX, 2)]; % Indices for edges (first and last row/column)
    for i = edgeIdx(1, :) % Horizontal edge lines (top and bottom)
        for j = 1:size(mesh.IX, 2) - 1
            nodeStart = mesh.IX(i, j, e);
            nodeEnd = mesh.IX(i, j + 1, e);
            edgeLinesX = [edgeLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
            edgeLinesY = [edgeLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
        end
    end
    for j = edgeIdx(2, :) % Vertical edge lines (left and right)
        for i = 1:size(mesh.IX, 1) - 1
            nodeStart = mesh.IX(i, j, e);
            nodeEnd = mesh.IX(i + 1, j, e);
            edgeLinesX = [edgeLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
            edgeLinesY = [edgeLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
        end
    end
end


axis equal;

text(1.06, 1.06, '(1, 1) ', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k');
text(-1.05, -1.05, '(-1, -1)', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', 'k');
text(1.1, 0, '\it{\xi}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k');
text(0.05, 1.1, '\it{\eta}  ', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k');

zero1 = [0.158, 0.516];
zero2 = [0.518, 0.156];
arrowlength = 0.75;
annotation('arrow', [zero1(1), zero1(1) + arrowlength], [zero1(2), zero1(2)], 'Color', 'k');
annotation('arrow', [zero2(1), zero2(1)], [zero2(2), zero2(2) + arrowlength], 'Color', 'k');
box off;
axis off;

legend('show', 'Location', 'southeastoutside');
fig.Position = [319.4000 537.8000 1000 420];

[X_pts, Y_pts] = meshgrid(GL_points, GL_points);
if N<=2
    plot(X_pts(:), Y_pts(:), 'x', 'MarkerSize', 6, 'DisplayName', 'GL Quadrature Points', 'MarkerEdgeColor',defaultColors(2,:)); % Red circles for quadrature points
else
    plot(mesh.X(:, 2), mesh.X(:, 3), 'x', 'MarkerSize', 6, 'DisplayName', 'GLL Quadrature Points', 'MarkerEdgeColor',defaultColors(4,:)); % Red circles for quadrature points
end
if N==1
        plot(10,10, 'x', 'MarkerSize', 6, 'DisplayName', 'GLL Quadrature Points', 'MarkerEdgeColor',defaultColors(4,:)); % Red circles for quadrature points
end
enhance_plot(0, fntsize, linethk, mrkSize, 0);

plot(edgeLinesX, edgeLinesY, 'Color', 'k', 'LineStyle', '-', 'LineWidth', linethk-1, 'HandleVisibility', 'off'); % Edge lines
plot(mesh.X(:, 2), mesh.X(:, 3), 'ko', 'MarkerSize', mrkSize, 'LineWidth', linethk-1, 'DisplayName', 'Nodes'); % Black diamonds for nodes

if N==1
    exportgraphics(gca, ['Figures\quadratureLEGEND.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

fig.Position = [651.4000 419.4000 432 430.4000];
legend off
exportgraphics(gca, ['Figures\meshLegend_N',num2str(N),'.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');


end

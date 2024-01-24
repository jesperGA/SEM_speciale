clear 
close all
clc

% Add necessary paths
addpath('FEM')
addpath('MESH')
addpath('PLOT')

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

% This is usually done in two steps:
% the interval [0,LX]*[0,LY] is divided into NELX*NELY quadrangular elements
% (that's basically a quadrangular finite element mesh),
% then each element is sub-grided with (P+1)^2 Gauss-Lobatto-Legendre points 
% (GLL nodes) where P is the polynomial order of the spectral elements.
% Actually, in this example the mesh is so simple 
% that we do it in a single step without storing intermediate data.

%**** Set here the parameters of the square box domain and mesh : ****
LX=1;	% side-length
NELX = 2; % number of elements along each side
study.N = 1; % polynomial degree (inside each element, along each direction)
%********

LY=LX;
NELY = NELX;
dxe = LX/NELX;
dye = dxe;
NEL = NELX*NELY;
NGLL = study.N+1; % number of GLL nodes per element

% Generate the Spectral Element mesh
% The domain is partitioned into elements,
% each element contains a cartesian GLL subgrid
[mesh.IX,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);

% Example data: node coordinates and numbers
node_coordinates = [x, y];
node_numbers = 1:size(node_coordinates, 1);

% element_nodes = [mesh.IX(:, :, 3),mesh.IX(:, :, 4)];
% xx = node_coordinates(element_nodes, 1);
% yy = node_coordinates(element_nodes, 2);
% ratios = 1 + 1/4 * sinpi(xx(yy==1))./(1/2);
% yy(yy==yy(21)) = 1/2+(yy(yy==yy(21))-1/2).*ratios;    
% yy(yy==yy(16)) = 1/2+(yy(yy==yy(16))-1/2).*ratios;    
% yy(yy==yy(11)) = 1/2+(yy(yy==yy(11))-1/2).*ratios;
% yy(yy==yy(6)) = 1/2+(yy(yy==yy(6))-1/2).*ratios;
% node_coordinates(element_nodes, 2) = yy;


% Plot nodes
figure;
plot(node_coordinates(:, 1), node_coordinates(:, 2), '.');
hold on;

% Display node numbers
for i = 1:length(node_numbers)
    % text(node_coordinates(i, 1), node_coordinates(i, 2), num2str(node_numbers(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
title('Node Numbers Plot');

% Add grid lines and customize tick marks
grid on;
xticks(unique(sort(node_coordinates(:, 1))));
% yticks(unique(sort(node_coordinates(:, 2))));
xtickformat('%.2f');
ytickformat('%.2f');

xlim([-0.1 1.1])
ylim([-0.1 1.4])

% Plot elements as boxes
for e = 1:size(mesh.IX,3)
    element_nodes = mesh.IX(:, :, e);
    x_lower = min(node_coordinates(element_nodes, 1));
    x_upper = max(node_coordinates(element_nodes, 1));
    y_lower = min(node_coordinates(element_nodes, 2));
    y_upper = max(node_coordinates(element_nodes, 2));

    % Plot box
    % rectangle('Position', [x_lower, y_lower, x_upper - x_lower, y_upper - y_lower], 'EdgeColor', 'r');

    % Display element number at the center of the box
    text((x_lower + x_upper) / 2, (y_lower + y_upper) / 2, num2str((e)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'r','FontSize',20);

    outer_nodes = [element_nodes(1,1:end), element_nodes(2:end,end)', element_nodes(end,end-1:-1:1), element_nodes(end-1:-1:1,1)'];
    plot(node_coordinates(outer_nodes,1),node_coordinates(outer_nodes,2),'k')
end
mesh.X = [(1:length(x))' x y];
mesh.Material = [1.1 1.2 1.3 1.4];
mesh.bound = [1 0 0];
mesh.bound = [2 0 0];

% Hold off to stop adding to the current plot
hold off;


opt = Controller(mesh, study);


%##########################################################################
%                           ADJ & JGA                                     %
%                     asbjorn@dyre-jespersen.dk                           %
%                                                                         %
%                 SEM implementation 2d - NS equation                     %
%##########################################################################

clear;
close all;
clc;

% Add necessary paths
addpath('FEM');
addpath('MESH');
addpath('PLOT');

%-------------------------------------------------------------------------%
%                             Physical domain and meshing                 %
%-------------------------------------------------------------------------%
% Parameters of the square box domain and mesh
LX = 1;             % Side-length of the box
NELX = 2;           % Number of elements along each side
LY = LX;
NELY = NELX;

F = @(X, Y) sin(X) .* exp(-Y);

for i = 1:4
    study.N = i;        % Polynomial degree inside each element
    
    [mesh] = mesh2D(LX, LY, NELX, NELY, study.N);
    
    % % Plot mesh
    % plotMesh2D(mesh)
    
    %-------------------------------------------------------------------------%
    %                             Solve                                       %
    %-------------------------------------------------------------------------%
    opt = Controller(mesh, study);
    for j=1:4
        error(i,j) = norm(abs(opt.U(mesh.IX(:,:,j),1) - F(mesh.X(mesh.IX(:,:,j),2),mesh.X(mesh.IX(:,:,j),3))),'inf');
    end
end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%
scatter3(mesh.X(:, 2), mesh.X(:, 3), opt.U);
hold on
scatter3(mesh.X(:,2),mesh.X(:,3),F(mesh.X(:,2),mesh.X(:,3)),'*')
hold on;
scatter3(mesh.X(mesh.IX(:,:,1),2),mesh.X(mesh.IX(:,:,1),3),opt.U(mesh.IX(:,:,1),1))

% Plot 'stiffnes' matrix
figure()
imagesc(opt.A); % This function scales the colors
colorbar

sum(sum(abs(opt.A)))
sum(sum(abs(opt.B)))
   % Plot mesh
    plotMesh2D(mesh)
figure()
semilogy(error)


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
LY = 1;
NELY = NELX;

F = @(X, Y) sin(X) .* exp(-Y);

for i = 2:5
    study.N = i;        % Polynomial degree inside each element
    
    [mesh] = mesh2D(LX, LY, NELX, NELY, study.N);

    % Plot mesh
    plotMesh2D(mesh)
    % plotMesh2Droenquistex(mesh)

    [mesh] = boundaryConditions(mesh, NELX, NELY);
    
    %-------------------------------------------------------------------------%
    %                             Solve                                       %
    %-------------------------------------------------------------------------%
    opt = Controller(mesh, study);
    error(i) = norm(abs(opt.U(:,1) - F(mesh.X(:,2),mesh.X(:,3))),'inf');

    % for j=1:4
    %     error(i,j) = norm(abs(opt.U(mesh.IX(:,:,j),1) - F(mesh.X(mesh.IX(:,:,j),2),mesh.X(mesh.IX(:,:,j),3))),'inf');
    % end
end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%
% plot nodal solution
plotNodalSolution(F,mesh.X,opt.U)

% % Plot 'stiffnes' matrix
% plotStiffnessMatrix(opt)
% 
% % Plot mesh
% plotMesh2D(mesh)
% 
% % Plot convergence
plotConvergence(error)

% Ne = 18;
% [meshEsol] = mesh2D(LX, LY, NELX, NELY, Ne);
% 
% opt.Xe = zeros(Ne+1, Ne+1, opt.nel);
% opt.Ye = zeros(Ne+1, Ne+1, opt.nel);
% opt.Ue = zeros(Ne+1, Ne+1, opt.nel);
% % Loop over elements and integrate
% for e = 1:opt.nel
%     X = mesh.X(mesh.IX(:,:,e),2);
%     Y = mesh.X(mesh.IX(:,:,e),3);
%     Xe = meshEsol.X(meshEsol.IX(:,:,e),2);
%     Ye = meshEsol.X(meshEsol.IX(:,:,e),3);
%     X = reshape(X, study.N+1, study.N+1);
%     Y = reshape(Y, study.N+1, study.N+1);
%     opt.Xe(:,:,e) = reshape(Xe, Ne+1, Ne+1);
%     opt.Ye(:,:,e) = reshape(Ye, Ne+1, Ne+1);
% 
%     [Ue(:,:,e)] = lagrange2D(X', Y', study.N, reshape(opt.U(mesh.IX(:,:,e),1),study.N+1,study.N+1)', opt.Xe(:,:,e)', opt.Ye(:,:,e)');
% 
%     % Ucheck = reshape(Ue(:,:,e),numel(Ue(:,:,e)),1)
%     % plotNodalSolution(F,meshEsol.X(meshEsol.IX(:,:,e),:),Ucheck)
% end




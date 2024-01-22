clear 
close all
clc

% Add necessary paths
addpath('FEM')
addpath('MESH')
addpath('PLOT')
addpath('TOPOPT')
addpath('Saved_values')

study.maxit=1; 
study.p = 1;

% Call the mesh routine 
L = pi;
ne = 2;

study.N = 5;

study.etype = 'SEM';

[mesh] = mesh1D(L, ne, study.N);

study.f = exp(mesh.X(:,2)) .* (cos(mesh.X(:,2))-sin(mesh.X(:,2)));

opt = Controller(mesh, study);



figure()
hold on
title(['ne = ',num2str(ne),', ', 'N = ',num2str(study.N),', '])
plot(mesh.X(:,2),opt.U,'o')
for e = 1:opt.nel
    plot(opt.elementX(e,:),opt.elementSolution(e,:),'Color','red')
end
% plot(linspace(mesh.X(1,2),mesh.X(end,2)),-sin(linspace(mesh.X(1,2),mesh.X(end,2))))
enhance_plot(0,0,0,0,0)
legend('Nodal solution','Element solution','','Analytic solution','location','best')
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
ne = 6;

study.N = 1;

study.etype = 'SEM';

[mesh] = mesh1D(L, ne, study.N);
mesh.X(:,2) = linspace(mesh.X(1,2),mesh.X(end,2),length(mesh.X))

study.f = exp(mesh.X(:,2)) .* (cos(mesh.X(:,2))-sin(mesh.X(:,2)));

opt = Controller(mesh, study);

study.N = ne/2;
ne=2;
[mesh2] = mesh1D(L, ne, study.N);
mesh.X = mesh2.X;
mesh.X(:,2) = linspace(mesh.X(1,2),mesh.X(end,2),length(mesh.X));
mesh.IX = mesh2.IX;
opt.nel = 2;

resolution = 1e5;
opt.elementSolution = zeros(opt.nel,resolution);
opt.elementX = zeros(opt.nel,resolution);

opt.U(1:end/2+1)=-1;
opt.U(end/2+1:end)=1;

% Loop over elements
for e=1:opt.nel

    % Get coordinates
    x = mesh.X(mesh.IX(e,2):mesh.IX(e,end-1),2);
    
    % Get element dofs
    edof=mesh.IX(e,2):mesh.IX(e,study.N+2);

    knots = mesh.X(edof,2)';

    ydata = opt.U(edof)';

    t = linspace(x(1),x(end),resolution);

    opt.elementX(e,:) = t;
    opt.elementSolution(e,:) = LagrangeFormInterpolation(knots,ydata,t);

end



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
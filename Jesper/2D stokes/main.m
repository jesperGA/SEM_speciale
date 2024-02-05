clear
close all
clc

%% addpath to FEA, MESH and VISUALIZATION (ParaView)
addpath('SEM')
addpath('MESH')
addpath('PLOT')

mat = [1.1,1.2,1.3,1.4;
    1.1,1.2,1.3,1.4];

% GLL = 2:1:13;
GLL = 5;
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
for order = 1:numel(GLL)
    n_GLL = GLL(order);

    L = pi;

    [xi,w,~] = lglnodes(n_GLL-1);
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;
    %% MESH
    [iglobV, xNV,yNV] = MeshBox_mod(2,2,2,2,n_GLL,1);
    [iglobP, xNP,yNP] = MeshBox_mod(2,2,2,2,n_GLL-2,2);

    mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP,yNP];
    mesh.IXv = iglobV;mesh.Xv = [(1:numel(xNV)).',xNV,yNV];
    %% Generate system matrices
    opt = [];
    [opt,study] = AssemblyQuad(mesh,opt,study);
    disp('Assembly done')
    %%
    % Modidy stiffness matrix for BCs
    k_org = opt.K;
    free = diag(opt.Null);
    opt.P = opt.P-k_org*opt.g;
    opt.P(~free) = opt.g(~free);
    opt.K = opt.Null'*opt.K*opt.Null - (opt.Null-speye(size(opt.Null)));
    opt.M = opt.Null'*opt.M*opt.Null - (opt.Null-speye(size(opt.Null)));
    % Solve static problem
    opt.U = opt.K \ (opt.P);
    %% Analytical solution
    % [Xn,Yn,Un] = meshgrid(xN,yN,opt.U);

    x = linspace(0,1);
    y = linspace(0,1.25);
    [X,Y] = meshgrid(x,y);
    xN = mesh.X(:,2);yN = mesh.X(:,3);
    sol = sin(X).*exp(-Y);
    sol_points = sin(xN).*exp(-yN);

    % [interp_data] = twoD_element_interpolator(mesh,mesh_interp, opt.U, xN, yN);
    % solution_plot(interp_data,n_interp)

    % figure()
    % % surf(X,Y,sol)
    % scatter3(xN,yN,opt.U,'r')
    % hold on
    % scatter3(xN,yN,sol_points,'*k')
    % legend('Numerical','Analytical')

    %% Error calc
    for i = 1:4
        point = mesh.IX(:,:,i);
        point = point(:);

        error(i,order) = norm(opt.U(point)-sol_points(point),inf);
    end

end


figure();semilogy(2*GLL+1,error,'-o','LineWidth',3)
grid on
%
% Setting the font size to 1
set(gca, 'FontSize', 18);
%
% % Adding labels with LaTeX interpreter
xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Element 1','Element 2','Element 3','Element 4','Interpreter','latex')

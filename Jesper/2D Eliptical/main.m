clear
close all
clc
mat = [1.1,1.2,1.3,1.4;
    1.1,1.2,1.3,1.4];

GLL = 2:1:11;
% GLL = 5;
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
for order = 1:numel(GLL)
    n_GLL = GLL(order);

    L = pi;

    [xi,w,~] = lglnodes(n_GLL-1);
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;
    %% MESH
    [iglob, xN,yN] = MeshBox(1,1,2,2,n_GLL,1);
    % [iglob2, xN2,yN2] = MeshBox(1,1,2,2,n_interp,2);
    %Create mesh like  Rønquist in Figure X. 4 elements in a 2x2 constallation.
    %and a sinus shaped top.
    mesh = modify_to_roenquist_mesh(xN,yN,iglob);
    % mesh_interp = modify_to_roenquist_mesh(xN2,yN2,iglob2);
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
    % for i = 1:4
    %     point = mesh.IX(:,:,i);
    %     point = point(:);

        error(order) = norm(opt.U-sol_points,inf);
    % end

end

df = readmatrix('roenquist_convergence_Poisson.csv','Delimiter',',');

figure();semilogy(2*GLL-1,error,'-o','LineWidth',3)
hold on
semilogy(df(:,1),df(:,2),'-*','LineWidth',3)
grid on
%
% Setting the font size to 1
set(gca, 'FontSize', 18);
%
% % Adding labels with LaTeX interpreter
xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);
legend('Numerical','Roenquist','Interpreter','latex')

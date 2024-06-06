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
addpath('MISC');


%-------------------------------------------------------------------------%
%                             Physical domain and meshing                 %
%-------------------------------------------------------------------------%
% Parameters of the square box domain and mesh
L = 2;
ne = 10;
study.N = 9;

mesh = mesh1D(L, ne, study.N);
mesh.X = mesh.X - 1;
% plotMesh1D(mesh);

% mesh.bound = [1]

opt = Controller(mesh, study); 
lambda1=2*sort(real(eig(full(opt.C),1i*full(opt.B))))/(pi*study.N*ne);
lambda1=lambda1(end/2:end)
lambda1 = lambda1([1 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 21 22 23 24 25 26 27 28 30 31 32 33 34 35 36 38 39 40 41 42 43 44 45 46 37 29 20 11 2])

ne = 90;
study.N = 1;

mesh = mesh1D(L, ne, study.N);
mesh.X = mesh.X - 1;
% plotMesh1D(mesh);

% mesh.bound = [1]

opt = Controller(mesh, study); 

lambda2=2*sort(real(eig(full(opt.C),1i*full(opt.B))))/(pi*study.N*ne);
lambda2=lambda2(end/2:end);
temp = lambda2;
temp(end/2+1:end)=lambda2(end:-2:2);
temp(1:end/2)=lambda2(1:2:end);
lambda2=temp;



% 
% % Set solver and example
% solvers = {'Direct','pcg','Uzawa','UzawamodJGA'};
% study.solver = solvers{2};
% examples = {'Roenquist','Roenquist_Poisson','Bercovier_1','Bercovier_2','LidDriven'};
% study.example = examples{1};
% study.unsteady = 1;
% if study.unsteady
% n_t = 50;
%     study.t = linspace(0,1,n_t);
% else
%     n_t = 2;
%     study.t = linspace(0,1,n_t);
% end
% 
% % Loop through polynomial orders
% polOrders = 7:7;
% 
% %-------------------------------------------------------------------------%
% %                             Physical domain and meshing                 %
% %-------------------------------------------------------------------------%
% % Parameters of the square box domain and mesh
% if strcmp(study.example,'Roenquist') || strcmp(study.example,'Roenquist_Poisson')
%     LX = 2; % Side-length of the box
% elseif strcmp(study.example,'Bercovier_1') || strcmp(study.example,'Bercovier_2') || strcmp(study.example,'LidDriven')
%     LX = 1;
% end
% NELX = 2; % Number of elements along each side
% LY = LX;
% NELY = NELX;
% 
% % Define functions based on example
% if strcmp(study.example,'Roenquist')
%     % Ronquist example
%     study.U1 = @(X1,X2) 1-X2.^2;
%     study.U2 = @(X1,X2) X2.*0;
%     study.P  = @(X1,X2) sinpi(X1).*sinpi(X2);
%     study.F1 = @(X1,X2) 2 + pi*cospi(X1).*sinpi(X2);
%     study.F2 = @(X1,X2)     pi*sinpi(X1).*cospi(X2);
% elseif strcmp(study.example,'Roenquist_Poisson')
%     % Ronquist Poisson example
%     study.U1 = @(X1,X2) sin(X1) .* exp(-X2);
%     study.U2 = @(X1,X2) X1 .* 0;
%     study.P  = @(X1,X2) X1 .* 0;
%     study.F1 = @(X1,X2) X1 .* 0;
%     study.F2 = @(X1,X2) X2 .* 0;
% elseif strcmp(study.example,'Bercovier_1')
%     % Bercovier (1)
%     study.U1 = @(X1,X2) -256 .* X1 .^2 .* (X1 - 1) .^ 2 .* X2 .* (X2 - 1) .*(2 .* X2 - 1);
%     study.U2 = @(X1,X2) study.U1(X2,X1); % Wrong sign!?
%     study.P  = @(X1,X2) X1 .* 0; 
%     study.F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
%                         2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2));
%     study.F2 = @(X1,X2) study.F1(X2,X1);
% elseif strcmp(study.example,'Bercovier_2')
%     % Bercovier (2)
%     study.U1 = @(X1,X2) -256 .* X1 .^2 .* (X1 - 1) .^ 2 .* X2 .* (X2 - 1) .*(2 .* X2 - 1);
%     study.U2 = @(X1,X2) -study.U1(X2,X1); 
%     study.P  = @(X1,X2) (X1-1/2) .* (X2-1/2);
%     study.F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
%                         2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2)) + X2 - 1/2;
%     study.F2 = @(X1,X2) study.F1(X2,X1);
% elseif strcmp(study.example,'LidDriven')
%     % Bercovier (2)
%     study.U1 = @(X1,X2) X2 == 1;
%     study.U2 = @(X1,X2) X1 .* 0;
%     study.P  = @(X1,X2) X1 .* 0;
%     study.F1 = @(X1,X2) X1 .* 0;
%     study.F2 = @(X1,X2) X2 .* 0;
% end
% 
% error_u = zeros(1, 13);
% error_p = zeros(1, 13);
% for i = polOrders
% 
%     study.N = i; % Polynomial degree inside each element
% 
%     % Generate mesh
%     [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);
% 
%     % Apply boundary conditions
%     [mesh] = boundaryConditions(study, mesh, NELX, NELY, study.U1, study.U2, study.P);
% 
%     %-------------------------------------------------------------------------%
%     %                             Solve                                       %
%     %-------------------------------------------------------------------------%
% 
%     opt = Controller(mesh, study);  
% 
% % maxp=max(opt.p);
% % temp=study.P(mesh.Xp(:,2),mesh.Xp(:,3));
% % maxtemp=max(temp);
% % dif=maxp-maxtemp;
% % opt.p=full(opt.p)-dif;
% % [full(opt.p) temp(:)];
% 
%     [error_u(i), error_p(i)] = calcError(mesh,study,opt);
% 
% end
% %-------------------------------------------------------------------------%
% %                             Plot solution                               %
% %-------------------------------------------------------------------------%
% 
% % fig = figure(1);
% % fig.Position = [1 49 2048 1.1568e+03];
% % % 
% % subplot(2,2,1)
% % Plot mesh
% % plotMesh2D(mesh)
% % if strcmp(study.example,'Roenquist')
%     plotMesh2Droenquistex(mesh)
% %     saveas(gcf,'Figures\2Dmesh','epsc')
% % elseif strcmp(study.example,'Bercovier_1')
% %     plotMesh2DBercovierex(mesh)
% %     legend off
% %     saveas(gcf,'Figures\2DmeshBercovier_1','epsc')
% % end
% plotNodalSolution(study.P,mesh.Xp,opt.p,'$p$')
% % legend('SEM','Analytical','Location','northeast')
% % saveas(gcf,'Figures\NodalSol_$p$','epsc')
% plotNodalSolution(study.U1,mesh.X,opt.u1,'$u_1$')
% plotNodalSolution(study.U2,mesh.X,opt.u2,'$u_2$')
% if strcmp(study.example,'Roenquist')
%     zlim([-1 1])
% end
% % legend('SEM','Analytical','Location','northeast')
% % saveas(gcf,['Figures\NodalSol_','$u_2$'],'epsc')
% 
% figure;
% for n = 1:n_t
% clf
% plotSol2D(mesh,opt.u1(:,n),opt.u2(:,n))
% pause(0.01)
% end
% 
% % plotPressure(mesh,study,opt)
% 
% % % Plot matrix
% % chegg(opt.LHS)
% 
% % resolution=1e2;
% % plotSolContour(mesh,study,opt,resolution,NELX, NELY)
% 
% if strcmp(study.example,'Roenquist')
    % Plot convergence
    load Roenquist_convection_op_eig_1D.csv
    load Roenquist_convection_op_eig_1D_FEM.csv
    fig = figure;
    fig.Position = [744 631.4000 923.4000 446.6000]
    plotC(mesh,lambda1,Roenquist_convection_op_eig_1D)
    plotC(mesh,lambda2,Roenquist_convection_op_eig_1D_FEM)
    plot([0 45],[0 1],'Color','k')
    legend('SEM $(N=9)$ - Our Results','SEM $(N=9)$ - R\o nquist ','Lin. FEM - Our Results','Lin. FEM - R\o nquist ','Exact Solution', ...
        'Interpreter', 'latex', 'FontSize', 18, 'Location','eastoutside' )
    enhance_plot(0, 20, 0, 8, 0);
    ylim([0 1.1999])
    xlim([0 45])
    xticks([0:5:45])
    saveas(gcf,'Figures\C_eigenvalues','epsc')
% elseif strcmp(study.example,'Bercovier_1')
%     fig = figure;
%     fig.Position = [744 636.2000 560 441.8000];
%         % Plot convergence
%     defaultColors = get(groot, 'DefaultAxesColorOrder');
%     semilogy(5:2:length(error_u)*2+2,error_u(2:end),'o-');
% 
%     % Labels and title
%     xlabel('$N_t$', 'Interpreter', 'latex', 'FontSize', 18);
%     ylabel('$\left\|u-u_h\right\|_{L^{\infty}, GL}$', 'Interpreter', 'latex', 'FontSize', 18);
% 
%     grid on;
%     % xlim([-0.1, 1.1]);
%     ylim([1e-16, 1]);
% 
%     % legend('$\left\|u-u_h\right\|_{L^{\infty}, GL}$ - Our Results','$\left\|u-u_h\right\|_{L^{\infty}, GL}$ - R\o nquist ','$\left\|p-p_h\right\|_{L^{\infty}, G}$ - Our Results','$\left\|p-p_h\right\|_{L^{\infty}, G}$ - R\o nquist ', ...
%     %     'Interpreter', 'latex', 'FontSize', 18, 'Location','southwest' )
%     enhance_plot(0, 20, 0, 0, 0);
%     saveas(gcf,'Figures\convergence_Bercovier_1','epsc')
% end

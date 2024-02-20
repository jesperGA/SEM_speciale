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

% Set solver and example
solvers = {'Direct','pcg','Uzawa','UzawamodJGA'};
study.solver = solvers{3};
examples = {'Roenquist','Roenquist_Poisson','Bercovier_1','Bercovier_2'};
study.example = examples{1};

% Loop through polynomial orders
polOrders = 3:3;

%-------------------------------------------------------------------------%
%                             Physical domain and meshing                 %
%-------------------------------------------------------------------------%
% Parameters of the square box domain and mesh
if strcmp(study.example,'Roenquist') || strcmp(study.example,'Roenquist_Poisson')
    LX = 2; % Side-length of the box
elseif strcmp(study.example,'Bercovier_1') || strcmp(study.example,'Bercovier_2')
    LX = 1;
end
NELX = 2; % Number of elements along each side
LY = LX;
NELY = NELX;

% Define functions based on example
if strcmp(study.example,'Roenquist')
    % Ronquist example
    study.U1 = @(X1,X2) 1-X2.^2;
    study.U2 = @(X1,X2) X2.*0;
    study.P  = @(X1,X2) sinpi(X1).*sinpi(X2);
    study.F1 = @(X1,X2) 2 + pi*cospi(X1).*sinpi(X2);
    study.F2 = @(X1,X2)     pi*sinpi(X1).*cospi(X2);
elseif strcmp(study.example,'Roenquist_Poisson')
    % Ronquist Poisson example
    study.U1 = @(X1,X2) sin(X1) .* exp(-X2);
    study.U2 = @(X1,X2) X1 .* 0;
    study.P  = @(X1,X2) X1 .* 0;
    study.F1 = @(X1,X2) X1 .* 0;
    study.F2 = @(X1,X2) X2 .* 0;
elseif strcmp(study.example,'Bercovier_1')
    % Bercovier (1)
    study.U1 = @(X1,X2) -256 .* X1 .^2 .* (X1 - 1) .^ 2 .* X2 .* (X2 - 1) .*(2 .* X2 - 1);
    study.U2 = @(X1,X2) study.U1(X2,X1); % Wrong sign!?
    study.P  = @(X1,X2) X1 .* 0; 
    study.F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
                        2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2));
    study.F2 = @(X1,X2) study.F1(X2,X1);
elseif strcmp(study.example,'Bercovier_2')
    % Bercovier (2)
    study.U1 = @(X1,X2) -256 .* X1 .^2 .* (X1 - 1) .^ 2 .* X2 .* (X2 - 1) .*(2 .* X2 - 1);
    study.U2 = @(X1,X2) -study.U1(X2,X1); 
    study.P  = @(X1,X2) (X1-1/2) .* (X2-1/2);
    study.F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
                        2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2)) + X2 - 1/2;
    study.F2 = @(X1,X2) study.F1(X2,X1);
end

error_u = zeros(1, 13);
error_p = zeros(1, 13);
for i = polOrders

    study.N = i; % Polynomial degree inside each element
    
    % Generate mesh
    [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);

    % Apply boundary conditions
    [mesh] = boundaryConditions(study, mesh, NELX, NELY, study.U1, study.U2, study.P);

    %-------------------------------------------------------------------------%
    %                             Solve                                       %
    %-------------------------------------------------------------------------%

    opt = Controller(mesh, study);  

% maxp=max(opt.p);
% temp=study.P(mesh.Xp(:,2),mesh.Xp(:,3));
% maxtemp=max(temp);
% dif=maxp-maxtemp;
% opt.p=full(opt.p)-dif;
% [full(opt.p) temp(:)];

    [error_u(i), error_p(i)] = calcError(mesh,study,opt);

end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%

% fig = figure(1);
% fig.Position = [1 49 2048 1.1568e+03];
% % 
% subplot(2,2,1)
% Plot mesh
% plotMesh2D(mesh)
if strcmp(study.example,'Roenquist')
    plotMesh2Droenquistex(mesh)
    saveas(gcf,'2DmeshN3','epsc')
elseif strcmp(study.example,'Bercovier_1')
    plotMesh2DBercovierex(mesh)
    legend off
    saveas(gcf,'2DmeshBercovier_1','epsc')
end
plotNodalSolution(study.P,mesh.Xp,opt.p,'$p$')
legend('SEM','Analytical','Location','northeast')
saveas(gcf,['NodalSol_$p1$'],'epsc')
plotNodalSolution(study.U1,mesh.X,opt.u1,'$u_1$')
plotNodalSolution(study.U2,mesh.X,opt.u2,'$u_2$')
if strcmp(study.example,'Roenquist')
    zlim([-1 1])
end
legend('SEM','Analytical','Location','northeast')
saveas(gcf,['NodalSol_','$u_2$'],'epsc')


% plotSol2D(mesh,opt)

% plotPressure(mesh,study,opt)

% % Plot matrix
% chegg(opt.LHS)

% resolution=1e2;
% plotSolContour(mesh,study,opt,resolution,NELX, NELY)

if strcmp(study.example,'Roenquist')
    % Plot convergence
    load Roenquist_p.csv
    load Roenquist_u.csv
    fig = figure;
    fig.Position = [517 573.8000 787 504.2000];
    plotConvergence(mesh,error_u,Roenquist_u)
    plotConvergence(mesh,error_p,Roenquist_p)
    legend('$\left\|u-u_h\right\|_{L^{\infty}, GL}$ - Our Results','$\left\|u-u_h\right\|_{L^{\infty}, GL}$ - R\o nquist ','$\left\|p-p_h\right\|_{L^{\infty}, G}$ - Our Results','$\left\|p-p_h\right\|_{L^{\infty}, G}$ - R\o nquist ', ...
        'Interpreter', 'latex', 'FontSize', 18, 'Location','southwest' )
    enhance_plot(0, 20, 0, 0, 0);
    saveas(gcf,'convergence','epsc')
elseif strcmp(study.example,'Bercovier_1')
    fig = figure;
    fig.Position = [744 636.2000 560 441.8000];
        % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    semilogy(5:2:length(error_u)*2+2,error_u(2:end),'o-');
    
    % Labels and title
    xlabel('$N_t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\left\|u-u_h\right\|_{L^{\infty}, GL}$', 'Interpreter', 'latex', 'FontSize', 18);

    grid on;
    % xlim([-0.1, 1.1]);
    ylim([1e-16, 1]);

    % legend('$\left\|u-u_h\right\|_{L^{\infty}, GL}$ - Our Results','$\left\|u-u_h\right\|_{L^{\infty}, GL}$ - R\o nquist ','$\left\|p-p_h\right\|_{L^{\infty}, G}$ - Our Results','$\left\|p-p_h\right\|_{L^{\infty}, G}$ - R\o nquist ', ...
    %     'Interpreter', 'latex', 'FontSize', 18, 'Location','southwest' )
    enhance_plot(0, 20, 0, 0, 0);
    saveas(gcf,'convergence_Bercovier_1','epsc')
end

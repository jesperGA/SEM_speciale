%##########################################################################
%                           ADJ & JGA                                     %
%                     asbjorn@dyre-jespersen.dk                           %
%                                                                         %
%                 SEM implementation 2d - NS equation                     %
%##########################################################################

clear;
close all;
clc;
warning('off','all')

% Add necessary paths
addpath('FEM');
addpath('MESH'); 
addpath('PLOT');
addpath('MISC');
addpath('nekmatlab');

% Set solver and example
solvers       = {'Direct','pcg','Uzawa','UzawamodJGA'};
timeStepping  = {'BDF1_AB3','BDF3_EX3'};
examples      = {'Roenquist','Roenquist_Poisson','Bercovier_1','Bercovier_2','LidDriven','Roenquist_NS','Pipe'};
study.timeInt = timeStepping{2};
study.solver  = solvers{1};
study.example = examples{7};

% convergence loop
convergence = {'spatial','temporal'};
study.convergence = convergence{1};

study.Re = 100;

if strcmp(study.convergence,'spatial')
    dt_vect = 1e-4;
elseif strcmp(study.convergence,'temporal')
    dt_vect = [1e-1 4e-2 2e-2 1e-2 5e-3 2e-3];
end

for j = 1:length(dt_vect)
    dt = dt_vect(j);
    if strcmp(study.convergence,'spatial')
        study.t = [0+dt:dt:0.1];
        % Loop through polynomial orders
        polOrders =5:5;
    elseif strcmp(study.convergence,'temporal')
        % Loop through time step sizes
        study.t = [0+dt:dt:1];
        polOrders =10;
    end
    
    %-------------------------------------------------------------------------%
    %                             Physical domain and meshing                 %
    %-------------------------------------------------------------------------%
    % Parameters of the square box domain and mesh
    if strcmp(study.example,'Roenquist') || strcmp(study.example,'Roenquist_Poisson') || strcmp(study.example,'Roenquist_NS') 
        LX = 2; % Side-length of the box
        NELX = 2; % Number of elements along each side
        LY = LX;
        NELY = NELX;
    elseif strcmp(study.example,'Bercovier_1') || strcmp(study.example,'Bercovier_2') || strcmp(study.example,'LidDriven')
        LX = 1;
        NELX = 2; % Number of elements along each side
        LY = LX;
        NELY = NELX;
    elseif  strcmp(study.example,'Pipe')
        LX = 5;
        LY = 1;
        NELX = 10; % Number of elements along each side
        NELY = 2;
    end

    
    % Define functions based on example
    study = createStudy(study);
    
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
    
        % if strcmp(study.convergence,'spatial')
        %     [error_u(i), error_p(i)] = calcError(mesh,study,opt);
        % elseif strcmp(study.convergence,'temporal')
        %     [error_u(j), error_p(j)] = calcError(mesh,study,opt);
        % end
    
    end
end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%
% 
% % plotMesh2D(mesh)
% plotMesh2Droenquistex(mesh)
% 

% %%%% plot PRESSURE
% figure()
% for n = 1:1width(opt.p)/50:width(opt.p);
%     plotNodalSolution(study.P,mesh.Xp,opt.p(:,n),'$p$',study,n)
% end
% % legend('SEM','Analytical','Location','southoutside')
% enhance_plot(0, 30, 100, 100, 0);
% % saveas(gcf,['Figures\nodalsol_p',num2str(n)],'epsc')
% % 


%%
%%%% plotSol2D VELOCITY
nodes=mesh.bound_u1(:,1);
figure;
for n = 1:round(length(study.t)/50):length(study.t)
    clf
    plotSol2D(mesh,opt.u1(:,n),opt.u2(:,n),study,study.t(n))
    title(sprintf('Time = %0.2f', study.t(n))); % Update title with the current time
end
% legend('','','','','SEM','Analytical','Location','southoutside','NumColumns',2)
% saveas(gcf,'Figures\plotSol2D','epsc')

% % resolution=1e2;
% % plotSolContour(mesh,study,opt,resolution,NELX, NELY)

% %%%% SAVE TO NEK
% fname = 'pipe';
% RefineTimes = 10;
% T = 1:10:length(study.t);
% save2nek5000(fname, T, opt,mesh,study,RefineTimes)
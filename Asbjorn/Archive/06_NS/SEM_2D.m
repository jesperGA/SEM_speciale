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

% Set solver and example
solvers = {'Direct','pcg','Uzawa','UzawamodJGA'};
study.solver = solvers{1};
timeStepping = {'BDF1_AB3','BDF3_EX3'};
study.timeInt = timeStepping{2};
examples = {'Roenquist','Roenquist_Poisson','Bercovier_1','Bercovier_2','LidDriven','Roenquist_NS','Pipe'};
study.example = examples{1};

% convergence loop
convergence = {'spatial','temporal'};
study.convergence = convergence{1};

study.unsteady = 1;
study.Re = 1;

if strcmp(study.convergence,'spatial')
    dt_vect = 1e-4;
elseif strcmp(study.convergence,'temporal')
    dt_vect = [1e-1 4e-2 2e-2 1e-2 5e-3 2e-3];
end

for j = 1:length(dt_vect)
    dt = dt_vect(j);
    if study.unsteady
        if strcmp(study.convergence,'spatial')
            study.t = [0+dt:dt:0.1];
            % Loop through polynomial orders
            polOrders =5:5;
        elseif strcmp(study.convergence,'temporal')
            % Loop through time step sizes
            study.t = [0+dt:dt:1];
            polOrders =10;
        end
    else
        study.t = [0 1];
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
    
plotMesh2D(mesh)

        % Apply boundary conditions
        [mesh] = boundaryConditions(study, mesh, NELX, NELY, study.U1, study.U2, study.P);
    
        %-------------------------------------------------------------------------%
        %                             Solve                                       %
        %-------------------------------------------------------------------------%
    
        opt = Controller(mesh, study);
    
        if strcmp(study.convergence,'spatial')
            [error_u(i), error_p(i)] = calcError(mesh,study,opt);
        elseif strcmp(study.convergence,'temporal')
            [error_u(j), error_p(j)] = calcError(mesh,study,opt);
        end
    
    end
end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%

% fig = figure(1);
% fig.Position = [1 49 2048 1.1568e+03];
% % 
% subplot(2,2,1)
% Plot mesh
plotMesh2D(mesh)
% if strcmp(study.example,'Roenquist')
    plotMesh2Droenquistex(mesh)
    % saveas(gcf,'Figures\2Dmesh_NS','epsc')
% elseif strcmp(study.example,'Bercovier_1')
%     plotMesh2DBercovierex(mesh)
%     legend off
%     saveas(gcf,'Figures\2DmeshBercovier_1','epsc')
% end

% plotNodalSolution(study.U2,mesh.X,opt.u2,'$u_1$',study)
if strcmp(study.example,'Roenquist')
    figure()
for n = 1:width(opt.p)/50:width(opt.p)
    plotNodalSolution(study.P,mesh.Xp,opt.p(:,n),'$p$',study,n)
        % title(sprintf('Time = %0.2f', study.t(n))); % Update title with the current time
    % legend('SEM','Analytical','Location','southoutside','NumColumns',2)
    % drawnow; % Update the plot
    % 
    % % Capture the plot as an image
    % frame = getframe(gcf);
    % img = frame2im(frame);
    % [imgInd, cmap] = rgb2ind(img, 256); % Convert the image to an indexed image
    % 
    % % Write to the GIF File
    % if n == 1
    %     imwrite(imgInd, cmap, 'timeseries_p.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    % else
    %     imwrite(imgInd, cmap, 'timeseries_p.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    % end
end
% legend('SEM','Analytical','Location','southoutside')
enhance_plot(0, 30, 100, 100, 0);
saveas(gcf,['Figures\nodalsol_p',num2str(n)],'epsc')
plotNodalSolution(study.U1,mesh.X,opt.u1,'$u_1$',study)
    zlim([-1 1])
end
legend('SEM','Analytical','Location','northeast')
saveas(gcf,['Figures\NodalSol_','$u_2$'],'epsc')
% 

%%%% MAKE GIF
nodes=mesh.bound_u1(:,1);
figure;
for n = 1:round(length(study.t)/50):length(study.t)
clf
plotSol2D(mesh,opt.u1(:,n),opt.u2(:,n),study,study.t(n))

    title(sprintf('Time = %0.2f', study.t(n))); % Update title with the current time
    % legend('','','','','SEM','Analytical','Location','southoutside','NumColumns',2)
    % drawnow; % Update the plot
    % 
    % % Capture the plot as an image
    % frame = getframe(gcf);
    % img = frame2im(frame);
    % [imgInd, cmap] = rgb2ind(img, 256); % Convert the image to an indexed image
    % 
    % % Write to the GIF File
    % if n == 1
    %     imwrite(imgInd, cmap, 'timeseries.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    % else
    %     imwrite(imgInd, cmap, 'timeseries.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    % end
end
legend('','','','','SEM','Analytical','Location','southoutside','NumColumns',2)
saveas(gcf,'Figures\plotSol2D','epsc')


% nodes=mesh.bound_u1(:,1);
% figure;
% for n = 1:round(length(study.t)/50):length(study.t)
% clf
% plotSol2D(mesh,opt.u1(:,n),opt.u2(:,n),study,study.t(n))
%     constant = 1.3;
%     xlim([-constant constant]);
%     ylim([-constant constant]);
% 
%     % title(sprintf('Time = %0.2f', study.t(n))); % Update title with the current time
%     drawnow; % Update the plot
% 
% end
% % legend('','','','','SEM','Analytical','Location','southoutside','NumColumns',2)
% saveas(gcf,['Figures\plotSol2D',num2str(n)],'epsc')



% plotPressure(mesh,study,opt)

% % Plot matrix
% chegg(opt.LHS)

% resolution=1e2;
% plotSolContour(mesh,study,opt,resolution,NELX, NELY)
defaultColors = get(groot, 'DefaultAxesColorOrder');
if strcmp(study.example,'Roenquist_NS') && strcmp(study.convergence,'spatial')
    % Plot convergence
    load Roenquist_p_N.csv
    load Roenquist_u_N.csv
    fig = figure;
    fig.Position =[517 461 787 617];
    semilogy(5:2:length(error_u)*2+2,error_u(2:end),'o-','Color',defaultColors(1,:));
    hold on;
    semilogy([7 9 11],Roenquist_u_N(:,2),'*-','Color',defaultColors(2,:));
    semilogy(5:2:length(error_p)*2+2,error_p(2:end),'o-','Color',defaultColors(5,:));
    semilogy(Roenquist_p_N(:,1),Roenquist_p_N(:,2),'*-','Color',defaultColors(4,:));
    % Labels
    xlabel('$N_t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Error', 'Interpreter', 'latex', 'FontSize', 18);
    xlim([7 13])
    xticks([7 9 11 13])
    ylim([0.99e-5 1])
    grid on;
    legend('$\left\|u-u_h\right\|_{1, GL}$ - Our Results','$\left\|u-u_h\right\|_{1, GL}$ - R\o nquist ','$\left\|p-p_h\right\|_{1, G}$ - Our Results','$\left\|p-p_h\right\|_{1, G}$ - R\o nquist ', ...
        'Interpreter', 'latex', 'FontSize', 18, 'Location','eastoutside' )
    enhance_plot(0, 20, 0, 0, 0);
    saveas(gcf,'Figures\convergence_spatial','epsc')
elseif strcmp(study.example,'Roenquist_NS') && strcmp(study.convergence,'temporal')
    % Plot convergence
    load Roenquist_p.csv
    load Roenquist_u.csv
    fig = figure;
    fig.Position = [517 537.8000 787 540.2000];
    loglog(dt_vect,error_u,'o-','Color',defaultColors(1,:));
    hold on
    loglog(Roenquist_u(:,1),Roenquist_u(:,2),'*-','Color',defaultColors(2,:));    
    loglog(dt_vect,error_p,'o-','Color',defaultColors(5,:));
    loglog(Roenquist_p(:,1),Roenquist_p(:,2),'*-','Color',defaultColors(4,:));
    grid on;
    % Labels
    xlabel('$\Delta t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Error', 'Interpreter', 'latex', 'FontSize', 18);
    legend('$\left\|u-u_h\right\|_{1, GL}$ - Our Results','$\left\|u-u_h\right\|_{1, GL}$ - R\o nquist ','$\left\|p-p_h\right\|_{1, G}$ - Our Results','$\left\|p-p_h\right\|_{1, G}$ - R\o nquist ', ...
        'Interpreter', 'latex', 'FontSize', 18, 'Location','eastoutside' )
    enhance_plot(0, 20, 0, 0, 0);
    saveas(gcf,'Figures\convergence_temporal','epsc')
end


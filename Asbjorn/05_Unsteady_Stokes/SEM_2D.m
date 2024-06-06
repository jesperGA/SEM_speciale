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
study.solver = solvers{2};
examples = {'Roenquist','Roenquist_Poisson','Bercovier_1','Bercovier_2','LidDriven'};
study.example = examples{1};
study.unsteady = 1;
if study.unsteady
    dt = 0.01;
    study.t = [0:dt:0.5];
else
    n_t = 2;
    study.t = linspace(0,1,n_t);
end

% Loop through polynomial orders
polOrders = 3:3;

%-------------------------------------------------------------------------%
%                             Physical domain and meshing                 %
%-------------------------------------------------------------------------%
% Parameters of the square box domain and mesh
if strcmp(study.example,'Roenquist') || strcmp(study.example,'Roenquist_Poisson')
    LX = 2; % Side-length of the box
elseif strcmp(study.example,'Bercovier_1') || strcmp(study.example,'Bercovier_2') || strcmp(study.example,'LidDriven')
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
elseif strcmp(study.example,'LidDriven')
    % Bercovier (2)
    study.U1 = @(X1,X2) X2 == 1;
    study.U2 = @(X1,X2) X1 .* 0;
    study.P  = @(X1,X2) X1 .* 0;
    study.F1 = @(X1,X2) X1 .* 0;
    study.F2 = @(X1,X2) X2 .* 0;
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
% if strcmp(study.example,'Roenquist')
    plotMesh2Droenquistex(mesh)
%     saveas(gcf,'Figures\2Dmesh','epsc')
% elseif strcmp(study.example,'Bercovier_1')
%     plotMesh2DBercovierex(mesh)
%     legend off
%     saveas(gcf,'Figures\2DmeshBercovier_1','epsc')
% end


    figure()
for n = 1:width(opt.p)
    plotNodalSolution(study.P,mesh.Xp,opt.p(:,n),'$p$',study,n)
    %     title(sprintf('Time = %0.2f', study.t(n))); % Update title with the current time
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

    figure()
for n = 1:width(opt.u1)
    plotNodalSolution(study.U1,mesh.X,opt.u1(:,n),'$u_1$',study,n)
    %     title(sprintf('Time = %0.2f', study.t(n))); % Update title with the current time
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
    %     imwrite(imgInd, cmap, 'timeseries_u.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    % else
    %     imwrite(imgInd, cmap, 'timeseries_u.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    % end
end
% legend('SEM','Analytical','Location','southoutside')
enhance_plot(0, 30, 100, 100, 0);
saveas(gcf,['Figures\nodalsol_u',num2str(n)],'epsc')


% plotNodalSolution(study.U1,mesh.X,opt.u1,'$u_1$')
% plotNodalSolution(study.U2,mesh.X,opt.u2,'$u_2$')
% if strcmp(study.example,'Roenquist')
%     zlim([-1 1])
% end
% legend('SEM','Analytical','Location','northeast')
% saveas(gcf,['Figures\NodalSol_','$u_2$'],'epsc')

figure;
for n = 1:width(opt.u1);
clf
plotSol2D(mesh,opt.u1(:,n),opt.u2(:,n))
pause(0.01)
end

% plotPressure(mesh,study,opt)

% % Plot matrix
% chegg(opt.LHS)

% resolution=1e2;
% plotSolContour(mesh,study,opt,resolution,NELX, NELY)
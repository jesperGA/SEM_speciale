%##########################################################################
%                           ADJ & JGA                                     %
%                     asbjorn@dyre-jespersen.dk                           %
%                                                                         %
%                 SEM implementation 2d - NS equation                     %
%##########################################################################

clear;
close all;
% clc;
warning('off','all');

% Add necessary paths
addpath('FEM');
addpath('MESH'); 
addpath('PLOT');
addpath('MISC');

% nw=str2num(getenv('LSB_DJOB_NUMPROC'));
% pool = parpool(nw);
% gcp;

% Set solver and example configurations
timeSteppingMethods = {'BDF1_AB3','BDF3_EX3','FullyExplicit'};
exampleNames        = {'LidDriven','Roenquist_NS','Pipe'};
FEMTypes            = {'LFEM','QFEM','SEM'};

for femType=1:3
study.timeInt  = timeSteppingMethods{3};
study.example  = exampleNames{3};
study.FEM = FEMTypes{femType};
study.Re       = 100;
study.Save2NEK = 1;
% study.ObjectCoords = [0 0.25 0 0.5]; % [x1 x2 y1 y2] for square, [x, y, r] for circle
study.ObjectCoords = [0.2 0.2 0.1/2]; % [x1 x2 y1 y2] for square, [x, y, r] for circle
study.alpha = 1e4   ;
T = 8;


% Define convergence type and time step sizes
convergenceTypes = {'spatial','temporal'};
study.convergence = convergenceTypes{1};
dtValuesSpatial  =  1e-3;
dtValuesTemporal = [7e-3 4e-3 2e-3, 1e-3, 5e-4, 2e-4 1e-4];%, 5e-5, 2e-5, 1e-5];
dtValuesTemporal = [1e-3];%, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6];


dt_vect = dtValuesTemporal; % Default to temporal
if strcmp(study.convergence,'spatial')
    dt_vect = dtValuesSpatial;
end

i=0;
domain=2;
for dt = dt_vect

    study.t = 0 : dt : T;
    study.t_steps = 1:round(length(study.t)/1000):length(study.t);
    polOrders = 8; % Default for temporal
    if strcmp(study.convergence, 'spatial')
        polOrders = [12];
    end
    
    %-------------------------------------------------------------------------%
    %                             Physical domain and meshing                 %
    %-------------------------------------------------------------------------%
    
    % Define functions based on example
    study = createStudy(study);
    
    % Polynomial order loop
    for polOrder = polOrders

        study.N = polOrder; % Polynomial degree inside each element
        if strcmp(study.FEM,'LFEM')
            study.N = 1;
        elseif strcmp(study.FEM,'QFEM')
            study.N = 2;
        end

        % Setup physical domain and mesh
        [LX, LY, NELX, NELY] = defineDomain(study,polOrder,domain);

        % Generate mesh
        [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);

        plotMesh2Dgraphic(mesh,study.ObjectCoords)
legend off

        exportgraphics(gca, ['Figures\mesh',study.FEM,'polOrder',num2str(polOrder),'.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
        % plotMesh2D(mesh)

    end
end

end
        plotMesh2Dgraphic(mesh,study.ObjectCoords)
        exportgraphics(gca, ['Figures\meshLegend.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');


%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%

% plotMesh2D(mesh)
% plotMesh2Dgraphic(mesh,study.ObjectCoords)


% %%% plot PRESSURE
% figure()
% steps = 1:round(length(study.t)/100):length(study.t);
% for n = 1:width(opt.p)
%     time = study.t(steps(n));
%     plotNodalSolution(study.P,mesh.Xp,opt.p(:,n),'$p$',time)
% end
% enhance_plot(0, 30, 100, 100, 0);

% %%%% plot VELOCITY MAG
% figure()
% for n = 1:1width(opt.u1)/50:width(opt.p);
%     plotNodalSolution(study.P,mesh.Xp,opt.p(:,n),'$p$',study,n)
% end
% enhance_plot(0, 30, 100, 100, 0);

% %%
% %%%% plotSol2D VELOCITY
% nodes=mesh.bound_u1(:,1);
% fig=figure;
% fig.Position = [1 49 2048 1.1568e+03];
% for n = 1:10:find(sum(opt.u1),1,'last')
%     clf
%     time = study.t(study.t_steps(n));
%     plotSol2D(mesh,opt.u1(:,n),opt.u2(:,n),study,time)
%     title(sprintf('Time = %0.2f', time)); % Update title with the current time
% end


defaultColors = get(groot, 'DefaultAxesColorOrder');
if strcmp(study.example,'Roenquist_NS') && strcmp(study.convergence,'spatial')
    % Plot convergence
    load Roenquist_p_N.csv
    load Roenquist_u_N.csv
    fig = figure;
    fig.Position =[517 461 787 617];
    semilogy((polOrders+1)*2-1,error_u(error_u>0),'o-','Color',defaultColors(1,:));
    hold on;
    semilogy([7 9 11],Roenquist_u_N(:,2),'*-','Color',defaultColors(2,:));
    semilogy((polOrders+1)*2-1,error_p(error_p>0),'o-','Color',defaultColors(5,:));
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

function folderName = createSimulationFolder(mesh, study, dt, NELX, NELY)
    % createSimulationFolder Creates a folder for storing simulation results
    % Inputs:
    %   study - Struct with fields N, NELX, NELY that describe the simulation study
    %   dt    - Time step value as a double

    % Calculate dimensions based on study parameters
    NX = study.N * NELX + 1;
    NY = study.N * NELY + 1;

    % Determine the middle part of the folder name based on N
    if study.N == 1
        fem_type = 'LFEM';
    elseif study.N == 2
        fem_type = 'QFEM';
    else
        fem_type = 'SEM';
    end

    % Convert dt to a string in scientific notation and replace the decimal point
    dt_str = strrep(sprintf('%.1e', dt), '.', '_');

    % Create the full path for the new folder inside the parent folder
    if strcmp(study.convergence, 'spatial')
        folderName = fullfile('simulations2', sprintf('Vortex_%s_%dx%d_dt%s', fem_type, NX, NY, dt_str));
    elseif strcmp(study.convergence, 'temporal')
        folderName = fullfile('simulationsTemporal', sprintf('Vortex_%s_%dx%d_dt%s', fem_type, NX, NY, dt_str));
    end

    % Create the directory if it does not exist
    if ~exist(folderName, 'dir')
        mkdir(folderName);
        fprintf('Folder created: %s\n', folderName);
    end

    if study.Save2NEK
        x = mesh.X(mesh.X(:,3)==max(mesh.X(:,3)),2)';
        y = mesh.X(mesh.X(:,2)==max(mesh.X(:,2)),3)';
        filename = fullfile(folderName, 'meshNodeData.vtk');  % Output filename
        saveMeshToVTK(x, y, filename);
        
        cornerNodes = unique(mesh.IX([1 end], [1 end], :));
        x = unique(mesh.X(mesh.X(cornerNodes,1),2));
        y = unique(mesh.X(mesh.X(cornerNodes,1),3));
        filename = fullfile(folderName, 'meshElementData.vtk');  % Output filename
        saveMeshToVTK(x, y, filename);
    end


end



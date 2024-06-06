%##########################################################################
%                           ADJ & JGA                                     %
%                     asbjorn@dyre-jespersen.dk                           %
%                                                                         %
%                 SEM implementation 2d - NS equation                     %
%##########################################################################

% Clear workspace, close all figures, and clear command window
clear;
close all;
clc;
warning('off','all');

% Choose example problem
exampleNames = {'LidDriven', 'TaylorVortex', 'Duct'};
study.example = exampleNames{1}; % Choose example

% Add necessary paths
addpath('SEM', 'MESH', 'PLOT', 'MISC');

% Start parallel pool
try
    parpool('Threads', 8);
end

% Set solver method
timeSteppingMethods = {'BDF1_AB3', 'BDF3_AB3', 'FullyExplicit_FD1_AB3'}; % Semi-implicit methods
study.timeInt = timeSteppingMethods{3};

% Save simulation to NEK
study.Save2NEK = 1;

% Duct flow specific parameters
if strcmp(study.example, 'Duct')
    study.ObjectCoords = [0.2, 0.2, 0.1/2]; % [x, y, radius] for circle
    study.alpha = 0.1; % Brinkman penalty parameter
end

% Integration Time and Element polynomial order based on example
switch study.example
    case 'Duct'
        dt = 1e-4; % Time step size
        T = 8; % Total integration time
        study.N = 4; % Polynomial degree
    case 'LidDriven'
        dt = 1e-3; % Time step size
        T = 10; % Total integration time
        study.N = 10; % Polynomial degree
    case 'TaylorVortex'
        dt = 1e-3; % Time step size
        T = 0.5; % Total integration time
        study.N = 5; % Polynomial degree
end

study.t = 0:dt:T;
study.t_steps = 1:round(length(study.t)/1000):length(study.t); % Save roughly 1000 time steps

% Finite Element Method type
study.FEM = 'SEM'; % Default to Spectral Element Method
if study.N == 1
    study.FEM = 'LFEM';
elseif study.N == 2
    study.FEM = 'QFEM';
end

% Define functions based on example
study = createStudy(study);

% Setup physical domain and mesh
[LX, LY, NELX, NELY] = defineDomain(study);

% Generate mesh
mesh = meshStaggered(study, LX, LY, NELX, NELY, study.N);
if strcmp(study.example, 'Duct')
    [mesh.Object_nodes, mesh.scale] = insertObject(mesh, study);
end

% Plot mesh
plotMesh2Dgraphic(mesh, study);

% Apply boundary conditions
mesh = boundaryConditions(study, mesh, NELX, NELY, study.U1, study.U2, study.P);

%-------------------------------------------------------------------------%
%                    Assembly, Solve & Postprocess                        %
%-------------------------------------------------------------------------%

% Create simulation folder
study.folderName = createSimulationFolder(mesh, study, NELX, NELY)

% Initialize controller
opt = Controller(mesh, study);

%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%

if strcmp(study.example, 'TaylorVortex') || strcmp(study.example, 'LidDriven')
    % Plot velocity solution
    fig = figure;
    fig.Position = [1 49 2048 1.1568e+03];
    for n = 1:10:find(sum(opt.u1), 1, 'last')
        clf;
        time = study.t(study.t_steps(n));
        plotSol2D(mesh, opt.u1(:,n), opt.u2(:,n), study, time);
        title(sprintf('Time = %0.2f', time)); % Update title with the current time
    end
end

function folderName = createSimulationFolder(mesh, study, NELX, NELY)
    % Create a folder for storing simulation results

    % Create the full path for the new folder
    folderName = fullfile('simulations', sprintf('%s_%dx%d_Order%d', study.example, NELX, NELY, study.N));

    % Create the directory if it does not exist
    if ~exist(folderName, 'dir')
        mkdir(folderName);
        fprintf('Folder created: %s\n', folderName);
    end

    if study.Save2NEK
        % Save mesh to VTK format
        x = mesh.X(mesh.X(:,3) == max(mesh.X(:,3)), 2)';
        y = mesh.X(mesh.X(:,2) == max(mesh.X(:,2)), 3)';
        filename = fullfile(folderName, 'meshNodeData.vtk');
        saveMeshToVTK(x, y, filename);

        cornerNodes = unique(mesh.IX([1 end], [1 end], :));
        x = unique(mesh.X(mesh.X(cornerNodes, 1), 2));
        y = unique(mesh.X(mesh.X(cornerNodes, 1), 3));
        filename = fullfile(folderName, 'meshElementData.vtk');
        saveMeshToVTK(x, y, filename);
    end
end

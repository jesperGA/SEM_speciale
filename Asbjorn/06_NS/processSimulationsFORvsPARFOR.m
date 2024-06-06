%##########################################################################
%                           ADJ & JGA                                     %
%                     asbjorn@dyre-jespersen.dk                           %
%                                                                         %
%                 Load and analyze SEM simulations                        %
%##########################################################################

clear;
close all;
clc;

% Add necessary paths
addpath('FEM');
addpath('MESH'); 
addpath('PLOT');
addpath('MISC');

% Specify the root directory containing the folders
rootDir = 'c:\Users\asbjo\OneDrive - Danmarks Tekniske Universitet\Dokumenter\DTU\SEM Speciale\Kode\SEM_speciale\Asbjorn\06_NS\';

Dirr = [rootDir,'simulations_for\'];
allDataFor = loadData(Dirr)

Dirr = [rootDir,'simulations_parfor_Threads\'];
allDataThread = loadData(Dirr)
%%
CPUT_util_for = [17 40];
[LFEM_dofs_for, LFEM_t_assemblyC_for, SEM_dofs_for, SEM_t_assemblyC_for] = funktionen(allDataFor);

CPUT_util_parfor = [70 40];
[LFEM_dofs_parfor, LFEM_t_assemblyC_parfor, SEM_dofs_parfor, SEM_t_assemblyC_parfor] = funktionen(allDataThread);

% Plotting the combined data
plotTimeCombined(LFEM_dofs_for, LFEM_t_assemblyC_for, SEM_dofs_for, SEM_t_assemblyC_for, CPUT_util_for, ...
                 LFEM_dofs_parfor, LFEM_t_assemblyC_parfor, SEM_dofs_parfor, SEM_t_assemblyC_parfor, CPUT_util_parfor)

exportgraphics(gca, 'Figures\time_for_v_parfor.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

function plotTimeCombined(LFEM_dofs_for, LFEM_t_assemblyC_for, SEM_dofs_for, SEM_t_assemblyC_for, CPUT_util_for, ...
                          LFEM_dofs_parfor, LFEM_t_assemblyC_parfor, SEM_dofs_parfor, SEM_t_assemblyC_parfor, CPUT_util_parfor)

fig = figure;
fig.Position = [320 645 696.2000 315];
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order
colors = get(gca, 'ColorOrder');
fntSize = 18;

x_limits = [min([LFEM_dofs_for LFEM_dofs_parfor]) max([LFEM_dofs_for LFEM_dofs_parfor])]; % Get current x-axis limits to place the experimental value across the figure
xlim(x_limits)
ylim([0 100])

% Plot non-NaN data for 'for' loop
semilogx(LFEM_dofs_for, LFEM_t_assemblyC_for, 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', ['LFEM (for) - CPU: ',num2str(CPUT_util_for(1)),' %']);
semilogx(SEM_dofs_for, SEM_t_assemblyC_for, '^-', 'Color', colors(4,:), 'MarkerFaceColor', colors(4,:), 'DisplayName', ['SEM (for) - CPU: ',num2str(CPUT_util_for(2)),' %']);

% Plot non-NaN data for 'parfor' loop
semilogx(LFEM_dofs_parfor, LFEM_t_assemblyC_parfor, 'o--', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', ['LFEM (parfor) - CPU: ',num2str(CPUT_util_parfor(1)),' %']);
semilogx(SEM_dofs_parfor, SEM_t_assemblyC_parfor, '^--', 'Color', colors(4,:), 'MarkerFaceColor', colors(4,:), 'DisplayName', ['SEM (parfor) - CPU: ',num2str(CPUT_util_parfor(2)),' %']);

hold off;

% Add labels, legend, and grid
xlabel('Nr. of nodes in $x_1$-dir.', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;
enhance_plot(0, 20, 1.5, 8, 0);

end



function [LFEM_dofs, LFEM_t_assemblyC, SEM_dofs, SEM_t_assemblyC] = funktionen(data)

    % Assuming 'allData' is your main struct containing all the loaded data
    simulationNames = fieldnames(data);  % Get the list of all simulation names (keys in the struct)

    % Initialize containers for nodes and Strouhal numbers for each method
    LFEM_dofs = [];
    SEM_dofs = [];
    LFEM_t_assemblyC = [];
    SEM_t_assemblyC  = [];
    
    % Loop through each simulation
    for i = 1:length(simulationNames)
        simName = simulationNames{i};  % Get the simulation name
        currentOpt = data.(simName).opt;  % Access the 'opt' struct of the current simulation
        currentMesh = data.(simName).mesh;  % Access the 'opt' struct of the current simulation
        currentStudy = data.(simName).study;  % Access the 'opt' struct of the current simulation
    
        dofsX = sum(currentMesh.X(:,3)==min(currentMesh.X(:,3)));
    
        % Check the method type and store data accordingly
        if contains(simulationNames{i}, 'LFEM')
            LFEM_dofs = [LFEM_dofs, dofsX];
            LFEM_t_assemblyC = [LFEM_t_assemblyC, currentOpt.time_assemble_C];
        elseif contains(simulationNames{i}, 'SEM')
            SEM_dofs = [SEM_dofs, dofsX];
            SEM_t_assemblyC = [SEM_t_assemblyC, currentOpt.time_assemble_C];
        end
    end
end

function plotNodalTimeseries(X, U, time, nodeCoords, label, color)
    distances = sqrt((X(:,2) - nodeCoords(1)).^2 + (X(:,3) - nodeCoords(2)).^2);
    [~, loc] = min(distances);

    u = U(loc, :);
    hold on
    plot(time, u, 'DisplayName', label, 'Color', color,'LineWidth',1.5);
        findpeaks(u, time);
end

function allData = loadData(rootDir)

% Get a list of all subdirectories in the root directory
dirs = dir(rootDir);
dirs = dirs([dirs.isdir]);  % filter only directories
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));  % exclude current and parent directories

% Initialize a struct to hold all loaded data
allData = struct();

% Loop through each directory and load the .mat files
for i = 1:length(dirs)
    folderName = dirs(i).name;
    matFilePath = fullfile(rootDir, folderName, 'meshStudyOpt.mat');
    folderName = strrep(folderName,'-','m');

    % Check if the .mat file exists
    if exist(matFilePath, 'file')
        % Load the .mat file
        temp = load(matFilePath);
        
        % Check if required variables are in the file
        if isfield(temp, 'mesh') && isfield(temp, 'study') && isfield(temp, 'opt')
            % Assign to a new field in the allData struct
            allData.(folderName(8:end)).mesh = temp.mesh;
            allData.(folderName(8:end)).study = temp.study;
            allData.(folderName(8:end)).opt = temp.opt;
        else
            warning('File in %s does not contain the required variables.', folderName);
        end
    else
        warning('%s does not exist.', matFilePath);
    end
end

end

function St_values = calculateStrouhalPeaks(X, study, U, time, nodeCoords)
    D = 0.1; % Diameter
    U_m = study.Re / 100; % Mean velocity
    distances = sqrt((X(:,2) - nodeCoords(1)).^2 + (X(:,3) - nodeCoords(2)).^2);
    [~, loc] = min(distances);

    u = U(loc, :); % Velocity at the closest node

    % Find all peaks in the entire dataset
    [pks, locs] = findpeaks(u);
    pks = pks(2:2:end);
    locs=locs(2:2:end);
    if isempty(pks) || length(pks) < 2
        error('Not enough peaks found to perform analysis.');
    end

    % Initialize array to store Strouhal numbers
    St_values = zeros(1, length(pks) - 1); % Storage for up to the last 15 peaks, or fewer if not available

    % Loop over the number of last peaks to include from 1 to 15 (or max available)
    for numPeaks = 1:length(St_values)
        if numPeaks > 1 % Ensure at least two peaks to calculate a period
            last_pks_indices = locs(end-numPeaks+1:end); % Indices of the last 'numPeaks' peaks
            peak_times = time(last_pks_indices); % Corresponding times of these peaks
            periods = diff(peak_times); % Periods between successive peaks
            mean_period = mean(periods); % Average period
            frequency = 1 / mean_period; % Frequency calculation
            St = D * frequency / U_m; % Strouhal number calculation
        else
            St = NaN; % Cannot compute period with less than two peaks
        end

        St_values(numPeaks) = St; % Store computed Strouhal number
    end

    % Plot the Strouhal number as a function of number of peaks included
    figure;
    plot(1:length(St_values), St_values, '-o');
    xlabel('Number of Last Peaks Included');
    ylabel('Strouhal Number');
    title('Strouhal Number Calculation Based on Last Peaks Included');
    grid on;
end

function createAreaPlot(DoF, cMatrixTime, RHSTime, solvingTime, otherTime, isNaN, t_solveTot, fntSize, colors, figPosition, MrkSize)
    % Create the area plot
    fig = figure;
    fig.Position = figPosition;  % Set the position and size of the figure window
    hArea = area(DoF(~isNaN), [isNaN(~isNaN)' isNaN(~isNaN)' cMatrixTime(~isNaN)' isNaN(~isNaN)' RHSTime(~isNaN)' otherTime(~isNaN)' solvingTime(~isNaN)'], 'LineStyle', 'none');
    xlabel('Nr. of nodes in $x_1$-dir.', 'Interpreter', 'latex', 'FontSize', fntSize);
    ylabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', fntSize);
    xlim([min(DoF) max(DoF)]);
    ylim([0 15000]); % Set a fixed y-limit for consistency
    legend('','','Computing $\mathbf{C}^{n+1}\mathbf{u}^{n+1}$ $\phantom{quad}$','','Evaluating $\mathbf{g_i}\phantom{quad}$', 'Other$\phantom{quad}$', 'Solving $\mathbf{p}^{n+1}$ and $\mathbf{u}_i^{n+1}$', ...
        'Location', 'northoutside', 'Interpreter', 'latex', 'FontSize', fntSize,'NumColumns',4);
    grid on;
    hold on;

    % Plot NaN marks
    plot(DoF(isNaN), zeros(sum(isNaN), 1), 'x', 'MarkerSize', MrkSize, 'Color', colors(4,:), 'DisplayName', 'NaN - SEM ');

    % Plot total solve time
    plot(DoF(~isNaN), t_solveTot(~isNaN), '-o', 'Color', 'k', 'HandleVisibility', 'off');

    % Optional: Set the x-axis to logarithmic
    % set(gca, 'XScale', 'log');

    % Enhance plot appearance (adjust this function according to your settings)
    enhance_plot(0, 20, 1.5, MrkSize, 0);

    hold off;
end

function  [bw, type, E] = getBW(opt)
    H1 = opt.LHS_BC(1:opt.neqn_u,1:opt.neqn_u);
    H2 = opt.LHS_BC(opt.neqn_u+1:2*opt.neqn_u,opt.neqn_u+1:2*opt.neqn_u);
    D1 = -opt.LHS_BC(2*opt.neqn_u+1:end,1:opt.neqn_u);
    D2 = -opt.LHS_BC(2*opt.neqn_u+1:end,opt.neqn_u+1:2*opt.neqn_u);
    E = D1 * (H1 \ D1') + D2 * (H2 \ D2');
    bw = bandwidth(E);
        Dsysmat = decomposition(-E);
        type = Dsysmat.Type;
end

function [I, J, BD1, BD2] = precomputeGlobalCmatr(ldof, numElems, Be, D_hat, L1, L2, IX)
    % Precompute indices and coefficients for constructing global convection matrices.
    % The function computes tensor product derivatives D1 and D2, scales them with
    % element-specific length scales L1 and L2, and constructs sparse matrix indices.
    
    Imat = eye(length(D_hat));
    D1 = kron(Imat, D_hat) * 2;
    D2 = kron(D_hat, Imat) * 2;

    totalEntries = numElems * ldof^2;
    I = zeros(totalEntries, 1);
    J = zeros(totalEntries, 1);
    BD1 = zeros(totalEntries, 1);
    BD2 = zeros(totalEntries, 1);

    index = 0;
    for e = 1:numElems
        elementDofs = reshape(IX(:,:,e), [], 1);
        localIndices = (index + 1):(index + ldof^2);
        [rowIndices, colIndices] = ind2sub([ldof, ldof], localIndices - index);

        I(localIndices) = elementDofs(rowIndices);
        J(localIndices) = elementDofs(colIndices);

        % Compute scales for element-level convection matrix based on constant parts
        BD1_mat = Be(:,:,e) * (D1 / L1(e));
        BD2_mat = Be(:,:,e) * (D2 / L2(e));
        BD1(localIndices) = BD1_mat(:);
        BD2(localIndices) = BD2_mat(:);

        index = index + ldof^2;
    end
end

function C_temp = assembleC(I, J, ME_DE1_flat, ME_DE2_flat, U1, U2, neqnV)

    CE = U1(I) .* ME_DE1_flat + U2(I) .* ME_DE2_flat;
    
    % Assemble the sparse global matrix C
    C_temp = sparse(I, J, CE, neqnV, neqnV);
    
    % % Forming final output
    % Cv = [C_temp * U1;
    %       C_temp * U2;
    %       sparse(neqnP, 1)];
end
function [A, C, E] = getMatrices(mesh, study, opt)
% Initialize convection matrices
ldof = (study.N + 1)^2;
[~, ~, D_hat] = GetGLL(study.N+1);
[I, J, BD1, BD2] = precomputeGlobalCmatr(ldof, opt.nel, opt.b, D_hat, mesh.L1, mesh.L2, mesh.IX);
u1=opt.u1(:,end);
u2=opt.u2(:,end);
C = assembleC(I, J, BD1, BD2, full(u1), full(u2), opt.neqn_u);
[bw, type, E] = getBW(opt);
A = opt.A;
end

function plotSpy(A,C,E)

% Subplot 1 for matrix A
subplot(3, 1, 1);
spy(A);
title('$\mathbf{A}$', 'Interpreter', 'latex');
non_zeros_A = nnz(A);
non_zeros_A_text = sprintf('Non-zeros: %.0fK', non_zeros_A / 1000); % Format in thousands
text(0.5, -0.1, non_zeros_A_text, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
axis off; % Remove axis
enhance_plot(0, 0, 0, 0, 0);
% 
% % Subplot 2 for matrix C
% subplot(3, 1, 2);
% spy(C);
% title('$\mathbf{C}$', 'Interpreter', 'latex');
% non_zeros_C = nnz(C);
% non_zeros_C_text = sprintf('Non-zeros: %.0fK', non_zeros_C / 1000); % Format in thousands
% text(0.5, -0.1, non_zeros_C_text, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
% axis off; % Remove axis
% enhance_plot(0, 0, 0, 0, 0);

% Subplot 3 for matrix E
subplot(3, 1, 2);
spy(E);
title('$\mathbf{E}$', 'Interpreter', 'latex');
non_zeros_E = nnz(E);
non_zeros_E_text = sprintf('Non-zeros: %.0fK', non_zeros_E / 1000); % Format in thousands
text(0.5, -0.1, non_zeros_E_text, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
axis off; % Remove axis
enhance_plot(0, 0, 0, 0, 0);

% Adjust the figure
set(gcf, 'Position', [100, 100, 600, 800]); % Adjust figure size for better visualization
end

function [xDistance, yDistance] = findNearestXYDistances(points)
numPoints = size(points, 1);
xDistance = zeros(numPoints, 1);
yDistance = zeros(numPoints, 1);

% Extract unique x and y coordinates to determine grid spacing
uniqueX = unique(points(:,1));
uniqueY = unique(points(:,2));

% Calculate the minimum distances assuming uniform grid spacing
if length(uniqueX) > 1
    deltaX = min(diff(uniqueX)); % smallest x-axis spacing between neighbors
else
    deltaX = 0; % All points have the same x coordinate
end

if length(uniqueY) > 1
    deltaY = min(diff(uniqueY)); % smallest y-axis spacing between neighbors
else
    deltaY = 0; % All points have the same y coordinate
end

% Determine closest distances
for i = 1:numPoints
    currentPoint = points(i, :);

    % Find the closest point in the x direction
    if deltaX ~= 0
        xOptions = points(points(:,1) ~= currentPoint(1), :); % Exclude the current point
        [~, idxX] = min(abs(xOptions(:,1) - currentPoint(1)));
        xDistance(i) = abs(xOptions(idxX, 1) - currentPoint(1));
    end

    % Find the closest point in the y direction
    if deltaY ~= 0
        yOptions = points(points(:,2) ~= currentPoint(2), :); % Exclude the current point
        [~, idxY] = min(abs(yOptions(:,2) - currentPoint(2)));
        yDistance(i) = abs(yOptions(idxY, 2) - currentPoint(2));
    end
end
end

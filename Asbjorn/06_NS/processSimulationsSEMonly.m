%##########################################################################
%                           ADJ & JGA                                     %
%                     asbjorn@dyre-jespersen.dk                           %
%                                                                         %
%                 Load and analyze SEM simulations                        %
%##########################################################################

clear;
close all;
% clc;

% Add necessary paths
addpath('FEM');
addpath('MESH'); 
addpath('PLOT');
addpath('MISC');

% Specify the root directory containing the folders
rootDir = 'c:\Users\asbjo\OneDrive - Danmarks Tekniske Universitet\Dokumenter\DTU\SEM Speciale\Kode\SEM_speciale\Asbjorn\06_NS\SimulationsSEM\';

allData = loadData(rootDir)

St_bench = [0.2995	0.2941	0.2901	0.2979	0.296	0.2985	0.2959	0.2907	0.2967	0.333	0.338	0.3003	0.2973	0.2881	0.2994	0.2968	0.302	0.303	0.313	0.3	0.289	0.3017	0.3012	0.2957	0.2997	0.2927	0.2713	0.2768	0.2778	0.2646	0.2841	0.2336	0.274];

%%

% Assuming 'allData' is your main struct containing all the loaded data
simulationNames = fieldnames(allData);  % Get the list of all simulation names (keys in the struct)

numCols = length(simulationNames) / 4;

% Initialize the matrices
dofs       = zeros(4, numCols);
St         = zeros(4, numCols);
t_solveTot = zeros(4, numCols);
last_u1    = zeros(4, numCols);
N          = zeros(4, numCols);
K          = zeros(4, numCols);

% Loop through each simulation
for i = 1:length(simulationNames)
    simName = simulationNames{i};  % Get the simulation name
    currentOpt = allData.(simName).opt;  % Access the 'opt' struct of the current simulation
    currentMesh = allData.(simName).mesh;  % Access the 'opt' struct of the current simulation
    currentStudy = allData.(simName).study;  % Access the 'opt' struct of the current simulation

    if abs(0.29-currentOpt.St*2)<0.03
        currentOpt.St = currentOpt.St*2;
    elseif abs(0.29-currentOpt.St/3*2)<0.03
        currentOpt.St = currentOpt.St/3*2;
    end

    dofsX = sum(currentMesh.X(:,3)==min(currentMesh.X(:,3)));

    % [bw, type] = getBW(currentOpt);
    
    row = str2num(simName(end))-2;
    col = find(dofs(row,:)==0,1,'first');
    dofs(row,col) = dofsX;
    St(row,col) = currentOpt.St;
    t_solveTot(row,col) = currentOpt.solutionTime;
    last_u1(row,col) = currentOpt.u1(round(end/2),end);
    N(row,col) = currentStudy.N;
    K(row,col) = size(currentMesh.IX,3);
end

% Sort each row of St, t_solveTot, and last_u1 according to the corresponding row in dofs
for i = 1:size(dofs, 1)
    [dofs(i, :), sortedIndices] = sort(dofs(i, :));
    St(i, :) = St(i, sortedIndices);
    t_solveTot(i, :) = t_solveTot(i, sortedIndices);
    last_u1(i, :) = last_u1(i, sortedIndices);
    N(i, :) = N(i, sortedIndices);
end


fig = figure;
fig.Position = [319.4000 537.8000 405 360];
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order
colors = get(gca, 'ColorOrder');
fntSize = 18;
MrkSize = 12;
linethk = 2;
colorvec = colors([3 5 6 7],:);

% Add experimental Strouhal number
exp_St = 0.287;
exp_error = 0.003;
x_limits = [min(min(N)) max(max(N))]; % Get current x-axis limits to place the experimental value across the figure
% y_fill = [exp_St - exp_error, exp_St + exp_error];
% fill([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], [y_fill(1), y_fill(1), y_fill(2), y_fill(2)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Exp. St $\pm$ Error');

St_bench = [0.2995	0.2941	0.2901	0.2979	0.296	0.2985	0.2959	0.2907	0.2967	0.333	0.338	0.3003	0.2973	0.2881	0.2994	0.2968	0.302	0.303	0.313	0.3	0.289	0.3017	0.3012	0.2957	0.2997	0.2927	0.2713	0.2768	0.2778	0.2646	0.2841	0.2336	0.274];
median_bench = median(St_bench);
min_bench = min(St_bench);
max_bench = max(St_bench);

xlim(x_limits)
ylim([0.2 max_bench+0.01])

% Identify NaN indices for each method
isNaN = last_u1==0;

symbol = {'o-' 's-' '^-' 'd-'};
mesh = {'Very coarse' 'Coarse' 'Medium' 'Fine'};

% Plot non-NaN data
for i = 1:4
    semilogx(N(i,~isNaN(i,:)), St(i,~isNaN(i,:)), 'o-', 'Color', colorvec(i,:), 'DisplayName', mesh{i});
end
for i = 1:4
    % Highlight NaN entries with a special marker
    semilogx(N(i,isNaN(i,:)), 0.2+zeros(1, sum(isNaN(i,:))) + (i-1)*0.005, 'x', 'MarkerSize', MrkSize, 'Color', colorvec(i,:), 'DisplayName', ['NaN - ',mesh{i}]);
end


% Plot Benchmarks statistics
plot(x_limits, [median_bench, median_bench], '-', 'Color', 'k', 'LineWidth', linethk, 'DisplayName', 'Benchmarks Median');
plot(x_limits, [min_bench, min_bench], '--', 'Color', 'k','LineWidth',linethk , 'HandleVisibility',  'off');
plot(x_limits, [max_bench, max_bench], '--', 'Color', 'k', 'DisplayName', 'Benchmarks Min/Max');

plot(x_limits, [exp_St, exp_St], '-.', 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Exp. St');

hold off;

% Add labels, legend, and grid
xlabel('Polynomial order $N$', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Strouhal Number $St$', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;
enhance_plot(0, 20, linethk, 0, 0);
exportgraphics(gca, 'Figures\time_Convergence_St_SEM_LEGEND.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
legend off
exportgraphics(gca, 'Figures\Convergence_St_SEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%
fig = figure;
fig.Position = [319.4000 537.8000 405 360];
hold on; % Hold on to plot multiple datasets on the same figure

% % Get the default color order
% colors = get(gca, 'ColorOrder');
% fntSize = 18;

xlim(x_limits)
ylim([-max(max(t_solveTot))*0.09 max(max(t_solveTot))])

symbol = {'o-' 's-' '^-' 'd-'};
mesh = {'Very coarse' 'Coarse' 'Medium' 'Fine'};

% Plot non-NaN data
for i = 1:4
    semilogx(N(i,~isNaN(i,:)), t_solveTot(i,~isNaN(i,:)), 'o-', 'Color', colorvec(i,:), 'DisplayName', mesh{i});
end
for i = 1:4
    % Highlight NaN entries with a special marker
    semilogx(N(i,isNaN(i,:)), zeros(1, sum(isNaN(i,:))) + (i-4)*max(max(t_solveTot))*0.03, 'x', 'MarkerSize', MrkSize, 'Color', colorvec(i,:), 'DisplayName', ['NaN - ',mesh{i}]);
end
hold off;

% Add labels, legend, and grid
xlabel('Polynomial order $N$', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;
enhance_plot(0, 20, linethk, 8, 0);
legend off
exportgraphics(gca, 'Figures\time_Convergence_St_SEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%

% simName = 'SEM_82x55_dt1_0em04domain3';
% opt=allData.(simName).opt;
% study=allData.(simName).study;
% mesh=allData.(simName).mesh;
% 
% time = study.t(study.t_steps);
% vel_mag = sqrt(opt.u1.^2 + opt.u2.^2);
% nodeCoords = [0.615, 0.205];
% 
% 
% figure;
% % Plotting velocity vs. time
% plotNodalTimeseries(mesh.X, vel_mag, time, nodeCoords);
% % Plotting pressure vs. time
% plotNodalTimeseries(mesh.Xp, opt.p, time, nodeCoords);
% 
% % St_values = calculateStrouhalPeaks(mesh.X, study, vel_mag, time, nodeCoords)
% % St_value= calculateStrouhal(mesh.X, study, vel_mag, time, nodeCoords)


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
function plotNodalTimeseries(X, U, time, nodeCoords)
    distances = sqrt((X(:,2) - nodeCoords(1)).^2 + (X(:,3) - nodeCoords(2)).^2);
    [~, loc] = min(distances);

    u = U(loc, :);
    hold on
    plot(time, u);
    title('Velocity vs. Time');
    xlabel('Time (s)');
    ylabel('Velocity or Pressure (m/s) or (Pa)');
    findpeaks(u, time);
end
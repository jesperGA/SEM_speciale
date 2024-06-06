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
rootDir = 'c:\Users\asbjo\OneDrive - Danmarks Tekniske Universitet\Dokumenter\DTU\SEM Speciale\Kode\SEM_speciale\Asbjorn\06_NS\simulations\';

allData = loadData(rootDir)
%%
rootDir = 'c:\Users\asbjo\OneDrive - Danmarks Tekniske Universitet\Dokumenter\DTU\SEM Speciale\Kode\SEM_speciale\Asbjorn\06_NS\simulationsTemporal\';
allDataTemporal = loadData(rootDir)

St_bench = [0.2995	0.2941	0.2901	0.2979	0.296	0.2985	0.2959	0.2907	0.2967	0.333	0.338	0.3003	0.2973	0.2881	0.2994	0.2968	0.302	0.303	0.313	0.3	0.289	0.3017	0.3012	0.2957	0.2997	0.2927	0.2713	0.2768	0.2778	0.2646	0.2841	0.2336	0.274];
%%

% Initialize containers for nodes and Strouhal numbers for each method
LFEM_dofs = [];
LFEM_st = [];
QFEM_dofs = [];
QFEM_st = [];
SEM_dofs = [];
SEM_st = [];

LFEM_bw = [];
QFEM_bw = [];
 SEM_bw = [];
LFEM_type = {};
QFEM_type = {};
 SEM_type = {};

LFEM_t_assembly  = [];
LFEM_t_assemblyC = [];
LFEM_t_RHS       = [];
LFEM_t_BC        = [];
LFEM_t_solve     = [];
LFEM_t_solveTot  = [];
QFEM_t_assembly  = [];
QFEM_t_assemblyC = [];
QFEM_t_RHS       = [];
QFEM_t_BC        = [];
QFEM_t_solve     = [];
QFEM_t_solveTot  = [];
SEM_t_assembly   = [];
SEM_t_assemblyC  = [];
SEM_t_RHS        = [];
SEM_t_BC         = [];
SEM_t_solve      = [];
SEM_t_solveTot   = [];
LFEM_CFL = [];
QFEM_CFL = [];
SEM_CFL = [];



% Assuming 'allData' is your main struct containing all the loaded data
simulationNames = fieldnames(allData);  % Get the list of all simulation names (keys in the struct)

% Initialize a vector or cell array to store Strouhal numbers
StrouhalNumbers = [];  % Use a vector if all Strouhal numbers are numerical and always present

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

    [bw, type] = getBW(currentOpt);

    dt = currentStudy.t(2)-currentStudy.t(1);
    %Prepare for CFL calculations:
    [dx, dy] = findNearestXYDistances(currentMesh.X(:,2:3));
    if any(sum([currentOpt.u1])==0)
                CFL_max=NaN;
    else
        CFL = (currentOpt.u1.*dt./dx+currentOpt.u2.*dt./dy);
        CFL_max = (max(CFL));
            %         figure()
            %         semilogy(CFL_max)
            % hold on
            % title(currentStudy.folderName(25:end))
        CFL_max = max(CFL_max);
    end
    


    % Check the method type and store data accordingly
    if contains(simulationNames{i}, 'LFEM')
        LFEM_dofs = [LFEM_dofs, dofsX];
        LFEM_st = [LFEM_st, currentOpt.St];
        LFEM_t_assembly = [LFEM_t_assembly, currentOpt.assemblyTime];
        LFEM_t_assemblyC = [LFEM_t_assemblyC, currentOpt.time_assemble_C];
        LFEM_t_RHS = [LFEM_t_RHS, currentOpt.time_RHS];
        LFEM_t_BC  = [LFEM_t_BC, currentOpt.time_BC];
        LFEM_t_solve    = [LFEM_t_solve   , currentOpt.time_solve_p];
        LFEM_t_solveTot = [LFEM_t_solveTot, currentOpt.solutionTime];
        LFEM_bw = [LFEM_bw, bw];
        LFEM_type = {LFEM_type{:}, type};
        LFEM_CFL = [LFEM_CFL, CFL_max];
    elseif contains(simulationNames{i}, 'QFEM')
        QFEM_dofs = [QFEM_dofs, dofsX];
        QFEM_st = [QFEM_st, currentOpt.St];
        QFEM_t_assembly = [QFEM_t_assembly, currentOpt.assemblyTime];
        QFEM_t_assemblyC = [QFEM_t_assemblyC, currentOpt.time_assemble_C];
        QFEM_t_RHS = [QFEM_t_RHS, currentOpt.time_RHS];
        QFEM_t_BC  = [QFEM_t_BC, currentOpt.time_BC];
        QFEM_t_solve    = [QFEM_t_solve   , currentOpt.time_solve_p];
        QFEM_t_solveTot = [QFEM_t_solveTot, currentOpt.solutionTime];
        QFEM_bw = [QFEM_bw, bw];
        QFEM_type = {QFEM_type{:}, type};
        QFEM_CFL = [QFEM_CFL, CFL_max];
    elseif contains(simulationNames{i}, 'SEM')
        SEM_dofs = [SEM_dofs, dofsX];
        SEM_st = [SEM_st, currentOpt.St];
        SEM_t_assembly = [SEM_t_assembly, currentOpt.assemblyTime];
        SEM_t_assemblyC = [SEM_t_assemblyC, currentOpt.time_assemble_C];
        SEM_t_RHS = [SEM_t_RHS, currentOpt.time_RHS];
        SEM_t_BC  = [SEM_t_BC, currentOpt.time_BC];
        SEM_t_solve    = [SEM_t_solve   , currentOpt.time_solve_p];
        SEM_t_solveTot = [SEM_t_solveTot, currentOpt.solutionTime];
        SEM_bw = [SEM_bw, bw];
        SEM_type = {SEM_type{:}, type};
        SEM_CFL = [SEM_CFL, CFL_max];
    end
    StrouhalNumbers = [StrouhalNumbers, currentOpt.St];

end
%%
% Create a plot

fig = figure;
fig.Position = [320 620 700 340];
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order
colors = get(gca, 'ColorOrder');
fntSize = 18;
MrkSize = 12;

% Add experimental Strouhal number
exp_St = 0.287;
exp_error = 0.003;
x_limits = [min(LFEM_dofs) max(LFEM_dofs)]; % Get current x-axis limits to place the experimental value across the figure
% y_fill = [exp_St - exp_error, exp_St + exp_error];
% fill([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], [y_fill(1), y_fill(1), y_fill(2), y_fill(2)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Exp. St $\pm$ Error');

St_bench = [0.2995	0.2941	0.2901	0.2979	0.296	0.2985	0.2959	0.2907	0.2967	0.333	0.338	0.3003	0.2973	0.2881	0.2994	0.2968	0.302	0.303	0.313	0.3	0.289	0.3017	0.3012	0.2957	0.2997	0.2927	0.2713	0.2768	0.2778	0.2646	0.2841	0.2336	0.274];
median_bench = median(St_bench);
min_bench = min(St_bench);
max_bench = max(St_bench);

xlim(x_limits)
ylim([0.2 max_bench+0.01])

% Identify NaN indices for each method
LFEM_isNaN = isnan(LFEM_st);
QFEM_isNaN = isnan(QFEM_st);
SEM_isNaN = isnan(SEM_st);

% Plot non-NaN data
semilogx(LFEM_dofs(~LFEM_isNaN), LFEM_st(~LFEM_isNaN), 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', 'LFEM');
semilogx(QFEM_dofs(~QFEM_isNaN), QFEM_st(~QFEM_isNaN), 's-', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'DisplayName', 'QFEM');
semilogx( SEM_dofs(~SEM_isNaN), SEM_st(~SEM_isNaN), '^-', 'Color', colors(4,:), 'MarkerFaceColor',     colors(4,:), 'DisplayName', 'SEM');

% Highlight NaN entries with a special marker
semilogx(LFEM_dofs(LFEM_isNaN), zeros(sum(LFEM_isNaN), 1) + 0.22, 'x', 'MarkerSize', MrkSize, 'Color', colors(1,:), 'DisplayName', 'NaN - LFEM');
semilogx(QFEM_dofs(QFEM_isNaN), zeros(sum(QFEM_isNaN), 1) + 0.21, 'x', 'MarkerSize', MrkSize, 'Color', colors(2,:), 'DisplayName', 'NaN - QFEM');
semilogx( SEM_dofs(SEM_isNaN), zeros(sum(SEM_isNaN), 1)+0.2, 'x', 'MarkerSize', MrkSize,    'Color', colors(4,:), 'DisplayName', 'NaN - SEM');

% Plot Benchmarks statistics
plot(x_limits, [median_bench, median_bench], '-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Benchmarks Median');
plot(x_limits, [min_bench, min_bench], '--', 'Color', 'k','LineWidth',1.5 , 'HandleVisibility',  'off');
plot(x_limits, [max_bench, max_bench], '--', 'Color', 'k', 'DisplayName', 'Benchmarks Min/Max');

plot(x_limits, [exp_St, exp_St], '-.', 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Exp. St');

hold off;

% Add labels, legend, and grid
xlabel('Nr. of nodes in $x_1$-dir.', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Strouhal Number $St$', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;
enhance_plot(0, 20, 1.5, 0, 0);
exportgraphics(gca, 'Figures\Convergence_St_spatial.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
%%

%
fig = figure;
fig.Position = [320 620 420 340];
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order
colors = get(gca, 'ColorOrder');
fntSize = 18;

x_limits = [min(LFEM_dofs) max(LFEM_dofs)]; % Get current x-axis limits to place the experimental value across the figure
xlim(x_limits)
% ylim([0 max(StrouhalNumbers)])

% Plot non-NaN data
semilogx(LFEM_dofs(~LFEM_isNaN), LFEM_t_solveTot(~LFEM_isNaN), 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', 'LFEM');
semilogx(QFEM_dofs(~QFEM_isNaN), QFEM_t_solveTot(~QFEM_isNaN), 's-', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'DisplayName', 'QFEM');
semilogx( SEM_dofs(~SEM_isNaN),   SEM_t_solveTot(~SEM_isNaN), '^-', 'Color', colors(4,:), 'MarkerFaceColor',     colors(4,:), 'DisplayName', 'SEM');

% Highlight NaN entries with a special marker
semilogx(LFEM_dofs(LFEM_isNaN), zeros(sum(LFEM_isNaN), 1) + 0.02, 'x', 'MarkerSize', MrkSize, 'Color', colors(1,:), 'DisplayName', 'NaN - LFEM');
semilogx(QFEM_dofs(QFEM_isNaN), zeros(sum(QFEM_isNaN), 1) + 0.01, 'x', 'MarkerSize', MrkSize, 'Color', colors(2,:), 'DisplayName', 'NaN - QFEM');
semilogx( SEM_dofs(SEM_isNaN), zeros(sum(SEM_isNaN), 1), 'x', 'MarkerSize', MrkSize,    'Color', colors(4,:), 'DisplayName',       'NaN - SEM');

hold off;

% Add labels, legend, and grid
xlabel('Nr. of nodes in $x_1$-dir.', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;
enhance_plot(0, 20, 1.5, 8, 0);
legend off
exportgraphics(gca, 'Figures\time_Convergence_St_spatial.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%

%
fig = figure;
fig.Position = [319.4000 650.6000 351.2000 307.2000];
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order
colors = get(gca, 'ColorOrder');
fntSize = 18;

x_limits = [min(LFEM_dofs) max(LFEM_dofs)]; % Get current x-axis limits to place the experimental value across the figure
xlim(x_limits)
% ylim([0 max(StrouhalNumbers)])

% Plot non-NaN dat
plot(LFEM_dofs, LFEM_bw, 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', 'LFEM');
plot(QFEM_dofs, QFEM_bw, 's-', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'DisplayName', 'QFEM');
plot( SEM_dofs,  SEM_bw, '^-', 'Color', colors(4,:), 'MarkerFaceColor',     colors(4,:), 'DisplayName', 'SEM');

% % Highlight NaN entries with a special marker
% plot(LFEM_dofs(LFEM_isNaN), zeros(sum(LFEM_isNaN), 1) + 0.02, 'x', 'MarkerSize', MrkSize, 'Color', colors(1,:), 'DisplayName', 'LFEM NaN');
% plot(QFEM_dofs(QFEM_isNaN), zeros(sum(QFEM_isNaN), 1) + 0.01, 'x', 'MarkerSize', MrkSize, 'Color', colors(2,:), 'DisplayName', 'QFEM NaN');
% plot( SEM_dofs(SEM_isNaN), zeros(sum(SEM_isNaN), 1), 'x', 'MarkerSize', MrkSize,    'Color', colors(4,:), 'DisplayName', 'SEM NaN');

hold off;

% Add labels, legend, and grid
xlabel('Nr. of nodes in $x_1$-dir.', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Bandwidth $bw(\mathbf{E})$', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;
enhance_plot(0, 20, 1.5, 8, 0);
exportgraphics(gca, 'Figures\bw_converg.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%

% Define font size and colors
fntSize = 20;
colors = get(gca, 'ColorOrder');

% Define figure position
figPosition = [680 578 450 300];  % Modify as needed

% Call the function with LFEM data
createAreaPlot(LFEM_dofs, ...
               LFEM_t_assemblyC, ...
               LFEM_t_RHS, ...
               LFEM_t_solve, ...
               LFEM_t_solveTot - ...
               LFEM_t_assemblyC - ...
               LFEM_t_solve - ...
               LFEM_t_RHS, ...
               LFEM_isNaN, ...
               LFEM_t_solveTot, fntSize, colors, figPosition+[0 0 60 0],MrkSize);
legend off
exportgraphics(gca, 'Figures\timeLFEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% Call the function with LFEM data
createAreaPlot(QFEM_dofs, ...
               QFEM_t_assemblyC, ...
               QFEM_t_RHS, ...
               QFEM_t_solve, ...
               QFEM_t_solveTot - ...
               QFEM_t_assemblyC - ...
               QFEM_t_solve - ...
               QFEM_t_RHS, ...
               QFEM_isNaN, ...
               QFEM_t_solveTot, fntSize, colors, figPosition,MrkSize);
fig = gcf;
fig.Position = [91.4000 578 1.5784e+03 300];
exportgraphics(gca, 'Figures\timeLegend.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

fig.Position = figPosition
legend off
ylabel('')
yticklabels({'','','',''})
exportgraphics(gca, 'Figures\timeQFEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% Call the function with LFEM data
createAreaPlot(SEM_dofs, ...
               SEM_t_assemblyC, ...
               SEM_t_RHS, ...
               SEM_t_solve, ...
               SEM_t_solveTot - ...
               SEM_t_assemblyC - ...
               SEM_t_solve - ...
               SEM_t_RHS, ...
               SEM_isNaN, ...
               SEM_t_solveTot, fntSize, colors, figPosition,MrkSize);
for i = 1:length(SEM_dofs)
    if ~SEM_isNaN(i)
        if rem(i,2)==1
            text_pos = SEM_t_solveTot(i) + i^2 * 0.002 * max(SEM_t_solveTot);
        else
            text_pos = SEM_t_solveTot(i) + 0.1 * max(SEM_t_solveTot);
        end
        text(SEM_dofs(i), text_pos, SEM_type{i}(1:2), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', fntSize-4, 'Color', 'k');
        line([SEM_dofs(i) SEM_dofs(i)], [SEM_t_solveTot(i) text_pos], 'Color', colors(7,:));
    end
    end
    enhance_plot(0, 20, 2, MrkSize, 0);


legend off
ylabel('')
yticklabels({'','','',''})
exportgraphics(gca, 'Figures\timeSEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


%%

% Initialize containers for nodes and Strouhal numbers for each method
LFEM_dt = [];
LFEM_st = [];
LFEM_CFL = [];
QFEM_dt = [];
QFEM_st = [];
QFEM_CFL = [];
SEM_dt = [];
SEM_st = [];
SEM_CFL = [];

% Assuming 'allData' is your main struct containing all the loaded data
simulationNames = fieldnames(allDataTemporal);  % Get the list of all simulation names (keys in the struct)

% Initialize a vector or cell array to store Strouhal numbers
StrouhalNumbers = [];  % Use a vector if all Strouhal numbers are numerical and always present

% Loop through each simulation
for i = 1:length(simulationNames)
    simName = simulationNames{i};  % Get the simulation name
    currentOpt = allDataTemporal.(simName).opt;  % Access the 'opt' struct of the current simulation
    currentMesh = allDataTemporal.(simName).mesh;  % Access the 'opt' struct of the current simulation
    currentStudy = allDataTemporal.(simName).study;  % Access the 'opt' struct of the current simulation
    
    if abs(0.29-currentOpt.St*2)<0.03
        currentOpt.St = currentOpt.St*2;
    elseif abs(0.29-currentOpt.St/3*2)<0.03
        currentOpt.St = currentOpt.St/3*2;
    end
    dt = currentStudy.t(2)-currentStudy.t(1);

    %Prepare for CFL calculations:
    [dx, dy] = findNearestXYDistances(currentMesh.X(:,2:3));
    if any(sum([currentOpt.u1])==0)

                CFL_max=NaN;

    else
        CFL = (currentOpt.u1.*dt./dx+currentOpt.u2.*dt./dy);
        CFL_max = (max(CFL));
            %         figure()
            %         semilogy(CFL_max)
            % hold on
            % title(currentStudy.folderName(25:end))
        CFL_max = max(CFL_max);

    end
    
    % Check the method type and store data accordingly
    if contains(simulationNames{i}, 'LFEM')
        LFEM_dt = [LFEM_dt, dt];
        LFEM_st = [LFEM_st, currentOpt.St];
        LFEM_CFL = [LFEM_CFL, CFL_max];
    elseif contains(simulationNames{i}, 'QFEM')
        QFEM_dt = [QFEM_dt, dt];
        QFEM_st = [QFEM_st, currentOpt.St];
        QFEM_CFL = [QFEM_CFL, CFL_max];
    elseif contains(simulationNames{i}, 'SEM')
        SEM_dt = [SEM_dt, dt];
        SEM_st = [SEM_st, currentOpt.St];
        SEM_CFL = [SEM_CFL, CFL_max];
    end
    StrouhalNumbers = [StrouhalNumbers, currentOpt.St];

end

% Sort each row of St, t_solveTot, and last_u1 according to the corresponding row in dofs

    [LFEM_dt, sortedIndices] = sort(LFEM_dt);
     LFEM_st=                       LFEM_st( sortedIndices);
     LFEM_CFL=                       LFEM_CFL( sortedIndices);
    [QFEM_dt, sortedIndices] = sort(QFEM_dt);
     QFEM_st=                       QFEM_st( sortedIndices);  
     QFEM_CFL=                       QFEM_CFL( sortedIndices);
    [SEM_dt, sortedIndices] = sort(SEM_dt);
     SEM_st=                       SEM_st( sortedIndices);
     SEM_CFL=                       SEM_CFL( sortedIndices);

% Create a plot
%
fig = figure;
fig.Position = [319.4000 537.8000 810 420];
semilogx([1 2], [1e6 2e6], 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'HandleVisibility', 'off');
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order and set font size
colors = get(gca, 'ColorOrder');
fntSize = 18;

% Add experimental Strouhal number with experimental error shaded
exp_St = 0.287;
exp_error = 0.003;
y_fill = [exp_St - exp_error, exp_St + exp_error];
x_limits = [min([LFEM_dt, QFEM_dt, SEM_dt]), max([LFEM_dt, QFEM_dt, SEM_dt])]; % Adjusted to cover all datasets
% fill([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], [y_fill(1), y_fill(1), y_fill(2), y_fill(2)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Exp. St $\pm$ Error');

xlim(x_limits)
ylim([0.2 max_bench+0.01])

% Identify NaN indices for each method
LFEM_isNaN = isnan(LFEM_st);
QFEM_isNaN = isnan(QFEM_st);
SEM_isNaN = isnan(SEM_st);

% Plot non-NaN data
semilogx(LFEM_dt(~LFEM_isNaN), LFEM_st(~LFEM_isNaN), 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', 'LFEM');
semilogx(QFEM_dt(~QFEM_isNaN), QFEM_st(~QFEM_isNaN), 's-', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'DisplayName', 'QFEM');
semilogx(SEM_dt(~SEM_isNaN), SEM_st(~SEM_isNaN), '^-', 'Color', colors(4,:), 'MarkerFaceColor',     colors(4,:), 'DisplayName', 'SEM');

% Highlight NaN entries with a special marker
semilogx(LFEM_dt(LFEM_isNaN), zeros(sum(LFEM_isNaN), 1) + 0.22, 'x', 'MarkerSize', MrkSize, 'Color', colors(1,:), 'DisplayName', 'NaN - LFEM');
semilogx(QFEM_dt(QFEM_isNaN), zeros(sum(QFEM_isNaN), 1) + 0.21, 'x', 'MarkerSize', MrkSize, 'Color', colors(2,:), 'DisplayName', 'NaN - QFEM');
semilogx(SEM_dt(SEM_isNaN), zeros(sum(SEM_isNaN), 1)+0.2, 'x', 'MarkerSize', MrkSize,    'Color', colors(4,:), 'DisplayName', 'NaN - SEM');

St_bench = [0.2995	0.2941	0.2901	0.2979	0.296	0.2985	0.2959	0.2907	0.2967	0.333	0.338	0.3003	0.2973	0.2881	0.2994	0.2968	0.302	0.303	0.313	0.3	0.289	0.3017	0.3012	0.2957	0.2997	0.2927	0.2713	0.2768	0.2778	0.2646	0.2841	0.2336	0.274];
median_bench = median(St_bench);
min_bench = min(St_bench);
max_bench = max(St_bench);
% Plot benchmark statistics
plot(x_limits, [median_bench, median_bench], '-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Benchmarks Median');
plot(x_limits, [min_bench, min_bench], '--', 'Color', 'k','LineWidth',1.5 , 'HandleVisibility',  'off');
plot(x_limits, [max_bench, max_bench], '--', 'Color', 'k', 'DisplayName', 'Benchmarks Min/Max');

plot(x_limits, [exp_St, exp_St], '-.', 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Exp. St');


hold off;

% Add labels, legend, and grid
xlabel('Time step $\Delta t$', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Strouhal Number $St$', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;

% Function to enhance plot appearance (if defined)
enhance_plot(0, 20, 1.5, MrkSize, 0);

% Optionally save the plot
exportgraphics(gca, 'Figures\convergence_St_temporal.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


% Create a plot
%
fig = figure;
fig.Position =[319.4000 564.2000 730 400];
semilogx([1 2], [1e6 2e6], 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'HandleVisibility', 'off');
hold on; % Hold on to plot multiple datasets on the same figure

% Get the default color order and set font size
colors = get(gca, 'ColorOrder');
fntSize = 18;


x_limits = [min([LFEM_dt, QFEM_dt, SEM_dt]), max([LFEM_dt, QFEM_dt, SEM_dt])]; % Adjusted to cover all datasets
% fill([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], [y_fill(1), y_fill(1), y_fill(2), y_fill(2)], [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Exp. St $\pm$ Error');

% Identify NaN indices for each method
LFEM_isNaN = isnan(LFEM_st);
QFEM_isNaN = isnan(QFEM_st);
SEM_isNaN =   isnan(SEM_st);

xlim(x_limits)
ylim([0 max(max([LFEM_CFL(~LFEM_isNaN) QFEM_CFL(~QFEM_isNaN) SEM_CFL(~SEM_isNaN)]))+0.01])
ylim([0 1])


% Plot non-NaN data
semilogx(LFEM_dt(~LFEM_isNaN), LFEM_CFL(~LFEM_isNaN), 'o-', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'DisplayName', 'LFEM');
semilogx(QFEM_dt(~QFEM_isNaN), QFEM_CFL(~QFEM_isNaN), 's-', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'DisplayName', 'QFEM');
semilogx(SEM_dt(~SEM_isNaN),    SEM_CFL(~SEM_isNaN), '^-', 'Color', colors(4,:), 'MarkerFaceColor',     colors(4,:), 'DisplayName', 'SEM');

% Highlight NaN entries with a special marker
semilogx(LFEM_dt(LFEM_isNaN), zeros(sum(LFEM_isNaN), 1) + 0.06, 'x', 'MarkerSize', MrkSize, 'Color', colors(1,:), 'DisplayName', 'NaN - LFEM');
semilogx(QFEM_dt(QFEM_isNaN), zeros(sum(QFEM_isNaN), 1) + 0.03, 'x', 'MarkerSize', MrkSize, 'Color', colors(2,:), 'DisplayName', 'NaN - QFEM');
semilogx(SEM_dt(SEM_isNaN),    zeros(sum(SEM_isNaN), 1) + 0.0, 'x', 'MarkerSize', MrkSize,    'Color', colors(4,:), 'DisplayName', 'NaN - SEM');

hold off;

% Add labels, legend, and grid
xlabel('Time step $\Delta t$', 'Interpreter', 'latex', 'FontSize', fntSize);
ylabel('Courant Number $C   $', 'Interpreter', 'latex', 'FontSize', fntSize);
legend('show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', fntSize);
grid on;

% Function to enhance plot appearance (if defined)
enhance_plot(0, 20, 1.5, MrkSize, 0);

% Optionally save the plot
exportgraphics(gca, 'Figures\convergence_Courant.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%
simName = 'LFEM_91x61_dt1_0em04'
opt=allData.(simName).opt;
study=allData.(simName).study;
mesh=allData.(simName).mesh;
[A, C, E] = getMatrices(mesh, study, opt);
fig=figure;
plotSpy(A,C,E)
exportgraphics(fig, 'Figures\spyLFEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


simName = 'QFEM_91x61_dt1_0em04'
opt=allData.(simName).opt;
study=allData.(simName).study;
mesh=allData.(simName).mesh;
[A, C, E] = getMatrices(mesh, study, opt);
fig=figure;
plotSpy(A,C,E)
exportgraphics(fig, 'Figures\spyQFEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


simName = 'SEM_91x61_dt1_0em04'
opt=allData.(simName).opt;
study=allData.(simName).study;
mesh=allData.(simName).mesh;
[A, C, E] = getMatrices(mesh, study, opt);
fig=figure;
plotSpy(A,C,E)
exportgraphics(fig, 'Figures\spySEM.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');


simName = 'SEM_67x45_dt1_0em04'
opt=allData.(simName).opt;
study=allData.(simName).study;
mesh=allData.(simName).mesh;

time = study.t(study.t_steps);
vel_mag = sqrt(opt.u1.^2 + opt.u2.^2);
nodeCoords = [0.615, 0.205];

%%
fig=figure;
fig.Position = [744 689.8000 1000 300];
% Plotting velocity vs. time on the left y-axis
yyaxis left
plotNodalTimeseries(mesh.X, vel_mag, time, nodeCoords, 'Velocity', colors(1,:));
ylim([0 2])
ylabel('Velocity Mag.', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
enhance_plot(0, 20, 1.5, 0, 0);

% Plotting pressure vs. time on the right y-axis
yyaxis right
plotNodalTimeseries(mesh.Xp, opt.p, time, nodeCoords, 'Pressure', colors(2,:));
ylabel('Pressure', 'Interpreter', 'latex');
ylim ([0 2])
% legend({'Velocity', 'Pressure'}, 'Interpreter', 'latex');
% enhance_plot(0, 20, 1.5, 0, 0);
exportgraphics(gca, 'Figures\NodalTimeseries.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

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

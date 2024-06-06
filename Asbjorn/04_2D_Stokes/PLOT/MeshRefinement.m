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
study.solver = solvers{3};
examples = {'Roenquist','Roenquist_Poisson','Bercovier_1','Bercovier_2'};
study.example = examples{1};

% Loop through polynomial orders
polOrders = [2 4 6];

%-------------------------------------------------------------------------%
%                             Physical domain and meshing                 %
%-------------------------------------------------------------------------%
% Parameters of the square box domain and mesh
if strcmp(study.example,'Roenquist') || strcmp(study.example,'Roenquist_Poisson')
    LX = 2; % Side-length of the box
elseif strcmp(study.example,'Bercovier_1') || strcmp(study.example,'Bercovier_2')
    LX = 1;
end
NELX = 1; % Number of elements along each side
LY = LX;
NELY = NELX;

for i = polOrders

    study.N = i; % Polynomial degree inside each element
    
    % Generate mesh
    [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);


plotMesh2DroenquistexWOpressure(mesh,study.N)
saveas(gcf,['Figures\mesh','N',num2str(study.N),'nelX',num2str(NELX),'NELY',num2str(NELY)],'svg')

end

for i = [1 2 4 6]
    
    NELX = i; % Number of elements along each side
    NELY = NELX;

    study.N = 1; % Polynomial degree inside each element
    
    % Generate mesh
    [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);


plotMesh2DroenquistexWOpressure(mesh,study.N)
saveas(gcf,['Figures\mesh','N',num2str(study.N),'nelX',num2str(NELX),'nelY',num2str(NELY)],'svg')

end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%
for i = [1 2 3]
    
    NELX = i; % Number of elements along each side
    NELY = NELX;

    study.N = 2; % Polynomial degree inside each element
    
    % Generate mesh
    [mesh] = meshStaggered(study, LX, LY, NELX, NELY, study.N);


plotMesh2DroenquistexWOpressure(mesh,study.N)
saveas(gcf,['Figures\mesh','N',num2str(study.N),'nelX',num2str(NELX),'nelY',num2str(NELY)],'svg')

end
%-------------------------------------------------------------------------%
%                             Plot solution                               %
%-------------------------------------------------------------------------%




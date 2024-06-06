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

%-------------------------------------------------------------------------%
%                             Physical domain and meshing                 %
%-------------------------------------------------------------------------%
% Parameters of the square box domain and mesh
LX = 1;             % Side-length of the box
NELX = 2;           % Number of elements along each side
LY = 1;
NELY = NELX;

study.N = 3;
[mesh] = mesh2D(LX, LY, NELX, NELY, study.N);
plotMesh2D_coordinates(mesh)




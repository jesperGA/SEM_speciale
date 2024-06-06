function [mesh] = boundaryConditions(study, mesh, NELX, NELY, U1, U2, P)

LX = max(mesh.X(:,2)) - min(mesh.X(:,2));             % Side-length of the box
LY = max(mesh.X(:,3)) - min(mesh.X(:,3));
N = size(mesh.IX,1)-1;

%-------------------------------------------------------------------------%
%                             Boundary Conditions                         %
%-------------------------------------------------------------------------%

mesh.bound_u1= addBCs(U1,mesh,NELX,NELY,N);
mesh.bound_u2 = addBCs(U2,mesh,NELX,NELY,N);
mesh.bound_p = [];
% mesh.bound_p = [46,0,study.P(mesh.X(46,2),mesh.X(46,3))];

if strcmp(study.example,'Roenquist_Poisson')
    % Ronquist Poisson example
    mesh.bound_u2 = zeros(length(mesh.X),3);
    mesh.bound_u2(:,1) = mesh.X(:,1);
    mesh.bound_p = zeros(length(mesh.Xp),3);
    mesh.bound_p(:,1) = mesh.Xp(:,1);
elseif strcmp(study.example,'Bercovier_1') %|| strcmp(study.example,'LidDriven')
    mesh.bound_p = zeros(length(mesh.Xp),3);
    mesh.bound_p(:,1) = mesh.Xp(:,1);
end

end

function bounds = addBCs(F,mesh,NELX,NELY,N)
    % South, East, North, West BCs
    idx = 0;
    bounds = zeros((4*(NELX*N+1)-4),3);
    % South BCs
    idx=NELX*N+1;
    nodes = mesh.X(mesh.X(:,3)==min(mesh.X(:,3)),1);
    % East BCs
    idx=idx + NELY*N;
    temp = mesh.X(mesh.X(:,2)==max(mesh.X(:,2)),1); temp=temp(2:end); nodes = [nodes; temp];
    % North BCs
    idx=idx + NELX*N;
    temp = mesh.X(mesh.X(:,3)==max(mesh.X(:,3)),1); nodes=[nodes; temp(end-1:-1:1)];
    % West BCs
    idx=idx + NELY*N;
    temp = mesh.X(mesh.X(:,2)==min(mesh.X(:,2)),1);nodes=[nodes; temp(end-1:-1:2)];
    bounds(1:idx-1,1) =nodes;
    bounds(1:idx-1,3) = F(mesh.X(nodes,2),mesh.X(nodes,3));
end
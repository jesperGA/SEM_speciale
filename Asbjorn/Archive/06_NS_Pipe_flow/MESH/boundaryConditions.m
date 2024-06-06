function [mesh] = boundaryConditions(study, mesh, NELX, NELY, U1, U2, P)

LX = max(mesh.X(:,2)) - min(mesh.X(:,2));             % Side-length of the box
LY = max(mesh.X(:,3)) - min(mesh.X(:,3));
N = size(mesh.IX,1)-1;

%-------------------------------------------------------------------------%
%                             Boundary Conditions                         %
%-------------------------------------------------------------------------%

nodes = getBoundaryNodes(mesh);
if strcmp(study.example,'Pipe') & isfield(mesh,'newBoundary')
    % Example usage
    nodes = union(nodes,mesh.newBoundary);
end

mesh.bound_u1 = zeros(length(nodes),3);
mesh.bound_u2 = zeros(length(nodes),3);

mesh.bound_u1(:,1) = nodes;
mesh.bound_u2(:,1) = nodes;

mesh.bound_u1(:,3) = U1(mesh.X(nodes,2),mesh.X(nodes,3),0);
mesh.bound_u2(:,3) = U2(mesh.X(nodes,2),mesh.X(nodes,3),0);

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

function nodes = getBoundaryNodes(mesh,study)
    % South
    nodes = mesh.X(mesh.X(:,3)==min(mesh.X(:,3)),1);
    % East
    % temp = mesh.X(mesh.X(:,2)==max(mesh.X(:,2)),1); nodes = [nodes; temp(2:end)];
    % North
    temp = mesh.X(mesh.X(:,3)==max(mesh.X(:,3)),1); nodes=[nodes; temp(end:-1:1)];
    % West
    temp = mesh.X(mesh.X(:,2)==min(mesh.X(:,2)),1);nodes=[nodes; temp(end-1:-1:2)];
end
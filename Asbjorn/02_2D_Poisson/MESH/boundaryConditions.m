function [mesh] = boundaryConditions(mesh, NELX, NELY)

LX = max(mesh.X(:,2)) - min(mesh.X(:,2));             % Side-length of the box
LY = max(mesh.X(:,3)) - min(mesh.X(:,3));
N = size(mesh.IX,1)-1;

%-------------------------------------------------------------------------%
%                             Boundary Conditions                         %
%-------------------------------------------------------------------------%

F = @(X, Y) sin(X) .* exp(-Y);

mesh.bound= addBCs(F,mesh,NELX,NELY,N);

end

function bounds = addBCs(F,mesh,NELX,NELY,N)
    % South, East, North, West BCs
    idx = 0;
    bounds = zeros((4*(NELX*N+1)-4),3);
    % South BCs
    idx=NELX*N+1;
    nodes = mesh.X(mesh.X(:,3)==min(mesh.X(:,2)),1);
    % East BCs
    idx=idx + NELY*N;
    temp = mesh.X(mesh.X(:,2)==max(mesh.X(:,2)),1); temp=temp(2:end); nodes = [nodes; temp];
    % North BCs
    idx=idx + NELX*N;
    temp = [mesh.IX(:,end,3);mesh.IX(2:end,end,4)]; nodes=[nodes; temp(end-1:-1:1)];
    % West BCs
    idx=idx + NELY*N;
    temp = mesh.X(mesh.X(:,2)==min(mesh.X(:,3)),1);nodes=[nodes; temp(end-1:-1:2)];
    bounds(1:idx-1,1) =nodes;
    bounds(1:idx-1,3) = F(mesh.X(nodes,2),mesh.X(nodes,3));
end
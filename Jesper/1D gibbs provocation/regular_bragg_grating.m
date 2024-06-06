function mesh = regular_bragg_grating (L,N,gll,mat,study)
%mesh = regular_bragg_grating (L,N,mat)
% Generates mesh file for a bragg grate of L length that contains N number
% of nodes. With material alternating for each element.
nn = (N-1)*(gll-1)+1;
X = zeros(nn,2); %Generate nodes
X(:,1) = 1:nn;
% gll_fac = (L/N)*(study.xi+1)/2;
%
ne = N-1;
N = gll-1;
% IX = zeros(ne,1+gll);
%
%
% X(1,:) = [1 0];
% for i = 1:ne
%     node_ind = (i-1)*(gll-1)+2:i*(gll-1)+1;
%     X(node_ind,:) = [node_ind.',X(node_ind(1)-1,2) + gll_fac(2:end)];
%     IX(i,1:gll+1) = [i (i-1)*(gll-1)+1:i*(gll-1)+1];
%     if mod(i,2) == 0
%         IX(i,2+gll) = 2;
%     else
%         IX(i,2+gll) = 1;
%     end
% end
% mesh.IX = IX;
% mesh.X = X;
xi = study.xi;
mesh.X(1,1:2) = [1 0];
dx_e=L/ne;
for e=1:ne
    index_nodes=(e-1)*N+2:e*N+1;
    x_old=mesh.X(index_nodes(1)-1,2);
    mesh.X(index_nodes,1) = [(e-1)*(N)+2:e*N+1];
    mesh.X(index_nodes,2) = x_old+(xi(2:end)+1)/2*dx_e;
    edof = [(e-1)*N+1 e*N+1];
    mesh.IX(e,1)=e;
    mesh.IX(e,2:N+2)=edof(1):edof(2);

    if mod(e,2) == 0
        mesh.IX(e,2+gll) = 2;
    else
        mesh.IX(e,2+gll) = 1;
    end
end

% n_inactive = ne*2/39;
% n_middle =round( (ne - 2 * n_inactive) / n_div);
% 
% % Assign material 1 to the inactive ends
% mesh.IX(1:n_inactive, end) = 1;
% mesh.IX(end-n_inactive+1:end, end) = 1;
% 
% % Assign alternating materials 1 and 2 in the middle
% for i = 1:n_div
%     start_idx = n_inactive + (i - 1) * n_middle + 1;
%     end_idx = n_inactive + i * n_middle;
%     if mod(i, 2) == 1
%         mesh.IX(start_idx:end_idx, end) = 2;
%     else
%         mesh.IX(start_idx:end_idx, end) = 1;
%     end
% end

% material A E nu rho
mesh.Material = mat;
mesh.bound = [N,1,0];
mesh.PointLoads = [1,1,1];
% mesh.inactive = [1:n_inactive, (ne+1)-n_inactive:ne];

end
function [opt,study] = AssemblyQuad(mesh,opt,study)

n_GLL = study.n_GLL; %Specify number of GLL points
N = n_GLL-1;
ldof =(n_GLL^2);

opt.nel = size(mesh.IX,3);
opt.neqn = size(mesh.X,1)*2;

I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
KE = zeros(opt.nel*ldof*ldof,1);
ME = zeros(opt.nel*ldof*ldof,1);

ntriplets = 0;


%Define GLL weights and points. Collected from study struct.
w = study.w;
xi = study.xi;

for e = 1:opt.nel

    nen = mesh.IX(:,:,e);
    nen = nen(:);

    % generate edof
    % for i = 1:length(nen)
    %     edof(2*i-1) = 2*nen(i)-1;
    %     edof(2*i-0) = 2*nen(i)-0;
    % end
    %Calculate dofs for each elements. 
    edof = zeros(ldof,1);
    indices = 1:length(nen);
    edof(2*indices-1) = 2*nen - 1;
    edof(2*indices) = 2*nen;
    
    %Find x and y for each node. 
    xy = mesh.X(nen,2:3);
    x = reshape(xy(:,1),N+1,N+1);
    y = reshape(xy(:,2),N+1,N+1);
    % matID = mesh.IX(e,end);


    [me,ke] = twoD_helmholtz(x,y,[],[],n_GLL,w,xi);

    disp('test')
%     for krow = 1:ldof
%         for kcol = 1:ldof
%             ntriplets = ntriplets+1;
%             I(ntriplets) = edof(krow);
%             J(ntriplets) = edof(kcol);
%             KE(ntriplets) = ke(krow,kcol);
%             % mass matrix
%             ME(ntriplets) = me(krow,kcol);
% 
%         end
% 
%     end
% 
% end
% ind = find(I>0);
% opt.K = sparse(I(ind),J(ind),KE(ind),opt.neqn,opt.neqn);
% opt.M = sparse(I(ind),J(ind),ME(ind),opt.neqn,opt.neqn);

% Initial conditions for transient problem.

end
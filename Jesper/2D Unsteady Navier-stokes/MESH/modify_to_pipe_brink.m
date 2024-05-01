function mesh = modify_to_pipe_brink(h,L,mesh,V0)
%mesh = modify_to_pipe(h,L,mesh,V0,element_to_remove)
%Removes elements between 2<=x<=2.5 and 0<=y<=0.5

x = mesh.Xv(:,2); y = mesh.Xv(:,3);

body_bc = mesh.Xv(:,2) >=2 & mesh.Xv(:,2) <= 2.5 & mesh.Xv(:,3) <= 0.5;
mesh.body_ind = mesh.Xv(body_bc,1);

%
% [Xv_new, IXv_new, newBC] = remove_element(mesh.Xv,mesh.IXv,element_to_remove);
% [Xp_new, IXp_new, ~] = remove_element(mesh.Xp,mesh.IXp,element_to_remove);



west_bound = mesh.Xv(x==0,1);
north_bound = mesh.Xv(y==h,1);
south_bound = mesh.Xv(y==0,1);
% east_bound = mesh.Xv(x==L,1); %No boundaries at x=L



%Boundaries on west FOR V1:

west_bound_mat = [west_bound(:),ones(numel(west_bound),1),repmat(V0,numel(west_bound),1)];
south_bound_mat = [south_bound(:),ones(numel(south_bound),1),zeros(numel(south_bound),1)];
north_bound_mat = [north_bound(:),ones(numel(north_bound),1),zeros(numel(north_bound),1)];
mesh.bound = [west_bound_mat;south_bound_mat;north_bound_mat];

%Boundaries on west for V2

west_bound_mat = [west_bound(:),2*ones(numel(west_bound),1),zeros(numel(west_bound),1)];
south_bound_mat = [south_bound(:),2*ones(numel(south_bound),1),zeros(numel(south_bound),1)];
north_bound_mat = [north_bound(:),2*ones(numel(north_bound),1),zeros(numel(north_bound),1)];

mesh.bound = [mesh.bound;west_bound_mat;south_bound_mat;north_bound_mat];
% figure()
% plotmesh(iglob,xN,yN)

end
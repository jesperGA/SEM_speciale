function mesh = cylinder_brink(h,p,D,mesh,V0,fully)
%mesh = cylinder_brink(h,L,mesh,V0,element_to_remove)

x = mesh.Xv(:,2); y = mesh.Xv(:,3);

dist_vec = mesh.Xv(:,2:3)-p;
dist_norm = sqrt(sum(dist_vec.^2,2));

body_bc = dist_norm<=D/2;
mesh.body_ind = mesh.Xv(body_bc,1);

a0 = D/2;
a1 = D/4;
slope = (0-1)/(a0-a1);

fac = ones(numel(mesh.body_ind),2);
fac(:,1) = mesh.body_ind;
% Calculate the factors for points where distance is between D/3 and D/2
fac(:,2) = 1 + slope * (dist_norm(mesh.body_ind) - a1);
fac(dist_norm(mesh.body_ind) <= a1, 2) = 1;  % Ensure that distances <= D/3 have factor 1

mesh.alpha_fac = fac;

%Scale

%
% [Xv_new, IXv_new, newBC] = remove_element(mesh.Xv,mesh.IXv,element_to_remove);
% [Xp_new, IXp_new, ~] = remove_element(mesh.Xp,mesh.IXp,element_to_remove);



west_bound = mesh.Xv(x==0,1);
north_bound = mesh.Xv(y==h,1);
south_bound = mesh.Xv(y==0,1);
% east_bound = mesh.Xv(x==L,1); %No boundaries at x=L



%Boundaries on west FOR V1:
if ~fully
    west_bound_mat = [west_bound(:),ones(numel(west_bound),1),repmat(V0,numel(west_bound),1)];
elseif fully
    umax = 3/2*V0;
    R = h/2; r = linspace(0,2*R,numel(west_bound));
    U = umax.*(1-((r-R)/R).^2);
    west_bound_mat = [west_bound(:),ones(numel(west_bound),1),U'];
end

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
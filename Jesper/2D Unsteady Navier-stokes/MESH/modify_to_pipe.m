function mesh = modify_to_pipe(h,L,mesh,V0,element_to_remove,fully)
%mesh = modify_to_pipe(h,L,mesh,V0,element_to_remove)
%Removes elements between 2<=x<=2.5 and 0<=y<=0.5
%Find elements:
temp = 0;
for i = 1:size(mesh.IXv,3)

    nen = mesh.IXv(:,:,i);
    xx = mesh.Xv(nen,2);
    yy = mesh.Xv(nen,3);

    xxm = mean(xx);
    yym = mean(yy);

    if (xxm>2 && xxm<2.5) && yym<0.5
        if ~temp
            temp = 1;
            [Xv_new, IXv_new, ~] = remove_element(mesh.Xv,mesh.IXv,i);
            [Xp_new, IXp_new, ~] = remove_element(mesh.Xp,mesh.IXp,i);
        else
            [Xv_new, IXv_new, ~] = remove_element(Xv_new,IXv_new,i-temp);
            [Xp_new, IXp_new, ~] = remove_element(Xp_new,IXp_new,i-temp);
            temp = temp + 1;
        end

    end
end

mesh.IXp = IXp_new;mesh.Xp = Xp_new;
mesh.IXv = IXv_new;mesh.Xv = Xv_new;
x = mesh.Xv(:,2); y = mesh.Xv(:,3);

top_bc = mesh.Xv(:,2) >=2 & mesh.Xv(:,2) <= 2.5 & mesh.Xv(:,3) == 0.5;
west_bc = mesh.Xv(:,2) == 2 & mesh.Xv(:,3) <= 0.5;
east_bc = mesh.Xv(:,2) == 2.5 & mesh.Xv(:,3) <= 0.5;
box_BC = top_bc | west_bc | east_bc;
newBC = mesh.Xv(box_BC,1);

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
indent_bound_mat = [newBC,ones(numel(newBC),1),zeros(numel(newBC),1)];
mesh.bound = [west_bound_mat;south_bound_mat;north_bound_mat;indent_bound_mat];

%Boundaries on west for V2

west_bound_mat = [west_bound(:),2*ones(numel(west_bound),1),zeros(numel(west_bound),1)];
south_bound_mat = [south_bound(:),2*ones(numel(south_bound),1),zeros(numel(south_bound),1)];
north_bound_mat = [north_bound(:),2*ones(numel(north_bound),1),zeros(numel(north_bound),1)];
indent_bound_mat = [newBC,2*ones(numel(newBC),1),zeros(numel(newBC),1)];

mesh.bound = [mesh.bound;west_bound_mat;south_bound_mat;north_bound_mat;indent_bound_mat];
% figure()
% plotmesh(iglob,xN,yN)

end
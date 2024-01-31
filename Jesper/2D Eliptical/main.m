clear
close all
clc
mat = [1.1,1.2,1.3,1.4;
    1.1,1.2,1.3,1.4];

GLL = 2:1:13;
% GLL = 5;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
for order = 1:numel(GLL)
    n_GLL = GLL(order);

    L = pi;

    [xi,w,~] = lglnodes(n_GLL-1);
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;
    %% MESH
    [iglob, xN,yN] = MeshBox(1,1,2,2,n_GLL);
    %Create mesh like  Rønquist in Figure X. 4 elements in a 2x2 constallation.
    %and a sinus shaped top.
    for e = 3:4
        glob_mod = iglob(:,:,e);
        for k = 1:size(glob_mod,1)
            if e==3
                x_spec = xN(glob_mod(k,end));
                yN(glob_mod(k,end)) = 1+1/4*sin(pi*x_spec);
                ratio(k) = (1/2+1/4*sin(pi*x_spec))/(1/2);
                yN(glob_mod(k,2:end-1)) = (yN(glob_mod(k,2:end-1))-1/2)*(ratio(k))+1/2;
            elseif e==4 && k==1
            else
                x_spec = xN(glob_mod(k,end));
                yN(glob_mod(k,end)) = 1+1/4*sin(pi*x_spec);
                % ratio = (1/2-1/4*x_spec)/(1/2);
                yN(glob_mod(k,2:end-1)) = (yN(glob_mod(k,2:end-1))-1/2)*(ratio(end-k+1))+1/2;
            end

        end

    end

    %Boundaries on west
    west_bound = iglob(1,:,[1,3]);
    west_bound_mat = [west_bound(:),ones(numel(west_bound),1),zeros(numel(west_bound),1)];

    south_bound = iglob(:,1,[1,2]);
    south_bound_mat = [south_bound(:),ones(numel(south_bound),1),sin(xN(south_bound(:))).*exp(-yN(south_bound(:)))];

    east_bound = iglob(end,:,[2,4]);
    east_bound_mat = [east_bound(:),ones(numel(east_bound),1),sin(xN(east_bound(:))).*exp(-yN(east_bound(:)))];

    north_bound = iglob(:,end,[3,4]);
    north_bound_mat = [north_bound(:),ones(numel(north_bound),1),sin(xN(north_bound(:))).*exp(-yN(north_bound(:)))];

    mesh.bound = [west_bound_mat;south_bound_mat;east_bound_mat;north_bound_mat];
    % figure()
    % plotmesh(iglob,xN,yN)
    mesh.IX = iglob;
    mesh.X = [ones(length(xN),1),xN,yN]; %Save grid to mesh-struct
    %% Generate system matrices
    opt = [];
    [opt,study] = AssemblyQuad(mesh,opt,study);
    disp('Assembly done')
    %%
    % Modidy stiffness matrix for BCs
    k_org = opt.K;
    free = diag(opt.Null);
    opt.P = opt.P-k_org*opt.g;
    opt.P(~free) = opt.g(~free);
    opt.K = opt.Null'*opt.K*opt.Null - (opt.Null-speye(size(opt.Null)));
    opt.M = opt.Null'*opt.M*opt.Null - (opt.Null-speye(size(opt.Null)));
    % Solve static problem
    opt.U = opt.K \ (opt.P);
    %% Analytical solution
    % [Xn,Yn,Un] = meshgrid(xN,yN,opt.U);

    x = linspace(0,1);
    y = linspace(0,1.25);
    [X,Y] = meshgrid(x,y);

    sol = sin(X).*exp(-Y);
    sol_points = sin(xN).*exp(-yN);

    % figure()
    % % surf(X,Y,sol)
    % scatter3(xN,yN,opt.U,'r')
    % hold on
    % scatter3(xN,yN,sol_points,'*k')
    % legend('Numerical','Analytical')

    %% Error calc
    for i = 1
        point = mesh.IX(:,:,i);
        point = point(:);

        error(i,order) = norm(opt.U(point)-sol_points(point),inf);
    end

end


figure();semilogy(2*GLL+1,error,'-ok','LineWidth',3)
grid on
%
% Setting the font size to 18
set(gca, 'FontSize', 18);
%
% % Adding labels with LaTeX interpreter
xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);

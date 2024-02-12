clear
close all
clc

%% addpath to FEA, MESH and
addpath('SEM')
addpath('MESH')
addpath('PLOT')

mat = [1.1,1.2,1.3,1.4;
    1.1,1.2,1.3,1.4];

GLL = 4:1:14;
GLL = 6;
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
for order = 1:numel(GLL)
    n_GLL = GLL(order);

    L = pi;

    [xi,w,~] = lglnodes(n_GLL-1);
    [zeta,wp] = lgwt(n_GLL-2,-1,1);
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-2;
    study.zeta = zeta;study.wp = wp;
    %% MESH
    [iglobV, xNV,yNV] = MeshBox_mod(2,2,2,2,n_GLL,1);
    % mesh = modify_to_bercovier(xNV,yNV,iglobV);
    mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
    [iglobP, xNP,yNP] = MeshBox_mod(2,2,2,2,n_GLL-2,2);

    mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP-1,yNP-1];
    % mesh.IXv = iglobV;mesh.Xv = [(1:numel(xNV)).',xNV,yNV];

    %% Generate system matrices
    [opt, study] = controller(mesh, study)


    %Def P/force vector:

 


    u1 = opt.U(1:opt.neqnV);u2 = opt.U(opt.neqnV+1:2*opt.neqnV);p = opt.U(opt.neqnV*2+1:end);

    %% Analytical solution

    % U1 = @(xx,yy) -256 .* xx .^ 2 .* (xx - 1) .^ 2 .* yy .* (yy - 1) .* (2 .* yy - 1);
    % U2 = @(xx,yy) -U1(yy,xx);
    U1 = @(xx,yy) 1-yy.^2;
    U2 = @(xx,yy) 0;
    U_Tot = sqrt(U1(xx,yy).^2+U2(xx,yy).^2);
    u_mag = sqrt(u1.^2+u2.^2);
    xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
    sol1 = U1(xx,yy);
    sol2 = U2(xx,yy);
    psol = sin(pi.*xxp).*sin(pi*yyp);
    figure()
    scatter3(xx,yy,u1,'*k')
    hold on
    scatter3(xx,yy,sol1,'or')
    figure()
    scatter3(xx,yy,opt.P(opt.neqnV+1:2*opt.neqnV),'*k')

    figure()
    scatter3(xx,yy,u2,'*k')
    hold on
    scatter3(xx,yy,sol2,'or')

    figure()
    scatter3(xxp,yyp,p,'*k')
    hold on
    scatter3(xxp,yyp,psol,'or')
    % scatter3(xxp,yyp,zeros(size(psol)),'or')

    figure()
    scatter3(xx,yy,u_mag,'*k')
    hold on
    scatter3(xx,yy,U_Tot,'or')
   



    % % figure()
    % % % surf(X,Y,sol)
    % % scatter3(xN,yN,opt.U,'r')
    % % hold on
    % % scatter3(xN,yN,sol_points,'*k')
    % % legend('Numerical','Analytical')
    %
    % %% Error calc
    % for i = 1:4
    %     point = mesh.IX(:,:,i);
    %     point = point(:);
    %
        error(order) = norm(u_mag(:)-U_Tot(:),inf);
    % end

end

% A = readmatrix("Roenquist_u.csv");
% figure();semilogy(2*GLL-1,error,'-ok','LineWidth',3)
% hold on
% semilogy(A(:,1),A(:,2),'*r')
% grid on
% 
% % Setting the font size to 1
% set(gca, 'FontSize', 18);
% 
% % Adding labels with LaTeX interpreter
% xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);
% % legend('Element 1','Element 2','Element 3','Element 4','Interpreter','latex')
% legend('SEM','Roenquist','Interpreter','latex')


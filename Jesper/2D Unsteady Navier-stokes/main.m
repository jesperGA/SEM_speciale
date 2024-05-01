% clear
close all
clc

warning('off', 'all');

%% addpaths
addpath('SEM')
addpath('MESH')
addpath('PLOT')
addpath('misc')
%%
% mat = [1.1,1.2,1.3,1.4;
%     1.1,1.2,1.3,1.4];
% study.p_type = 'roenquist2';
% study.p_type = 'bercover';
study.p_type = 'liddriven';
study.solve_type = 'direct'; %uzawa
if strcmp(study.solve_type,'direct') == 1

    % study.direct_type = 'backslash';
    study.direct_type = 'LU';

end
%%
% study.solve_type = 'uzawa';
study.study_type = 'unsteady';
% study.precon = 'mhat';
study.precon = 'P';
% study.study_t ype = 'steady';
if strcmp(study.study_type,'unsteady') == 1
    study.T = 0.1;
    % study.nt = 1e3;
    % study.t = linspace(0,study.T,study.nt);
    study.t = 0:1e-3:study.T;
    study.nt = length(study.t);
    study.dt = (study.t(2)-study.t(1));

    % study.int_type = 'BDFk'; %Equivalent of solving Unsteady stokes.
    study.int_type = 'BDF3EX3'; %First order bdf for linear terms. 3 order for nonlinear terms.
    study.RE = 1000;
    study.BDF_order = 1;

    study.U10 = 0;
    study.U20 = 0;

    study.BC_type = 'dynamic';
    % study.BC_type = 'static';
end
%%
% GLL = 3:1:10;
GLL = 14;
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
figure(101)
dbp = readmatrix("Roenquist_p_N.csv");
dbu = readmatrix("Roenquist_u_N.csv");
% plot(db(:,2),db(:,1),'dk','DisplayName','Benchmark')
hold on

for order = 1:numel(GLL)
    n_GLL = GLL(order);


    error(order) = 0;
    errorp(order) = 0;

    [xi,w,~] = lglnodes(n_GLL-1);
    [zeta,wp] = lgwt(n_GLL-2,-1,1);
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-2;
    study.zeta = zeta;study.wp = wp;
    %% MESH
    [iglobV, xNV,yNV] = MeshBox_mod(5,5,1,1,n_GLL,1);
    % mesh = modify_to_bercovier(xNV,yNV,iglobV);
    % mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
    % mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
    mesh = liddriven(xNV,yNV,iglobV);
    [iglobP, xNP,yNP] = MeshBox_mod(5,5,1,1,study.n_GL,2);

    mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP,yNP];
    mesh.pref_dof = 1;

    %% Generate system matrices
    [opt, study] = controller(mesh, study);
    %Def P/force vector:
    u1 = opt.U(1:opt.neqnV);u2 = opt.U(opt.neqnV+1:2*opt.neqnV);p = opt.U(opt.neqnV*2+1:end);

    %% Analytical solution, validation and plots
    xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
    x2_ind = find(xx==0.5);
    xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);
    % U1 = @(xx,yy) -256 .* xx .^ 2 .* (xx - 1) .^ 2 .* yy .* (yy - 1) .* (2 .* yy - 1);
    % U2 = @(xx,yy) U1(yy,xx);
    U1 = @(xx,yy) 1-yy.^2;
    U2 = @(xx,yy) 0*xx;
    U_Tot = sqrt(U1(xx,yy).^2+U2(xx,yy).^2);
    u_mag = sqrt(u1.^2+u2.^2);
    p_sol = @(x,y,t) -1/4.*(cos(2.*x)+cos(2.*y)).*exp(-4.*t);

    % pGx = reshape(xxp,[2*study.n_GL,2*study.n_GL]);pGy = reshape(yyp,[2*study.n_GL,2*study.n_GL]);
    sol1 = U1(xx,yy);
    sol2 = U2(xx,yy);
    psol = sin(pi.*xxp).*sin(pi*yyp);
    % psol = zeros(length(xxp),1);
    % psol = (xxp-1/2).*(yyp-1/2);

    % u_check = [u1
    %     u2
    %     p];
    % P_check_vec = opt.sys_mat*u_check-opt.P;

    for i =1:size(opt.U,2)
        % figure(100)
        % clf
        % %     scatter3(xx,yy,opt.U(1:opt.neqnV,i),'*k')
        % %     scatter3(xxp,yyp,opt.Pr(:,i),'*k')
        sol = opt.g_sys(xx,yy,study.t(i));

        % plotSol2D(mesh,opt.U(1:opt.neqnV,i),opt.U(opt.neqnV+1:end,i),0)
        % hold on
        % plotSol2D(mesh,sol(1:opt.neqnV),sol(opt.neqnV+1:end),1)

        % mag = sqrt(opt.U(1:opt.neqnV,i).^2+opt.U(opt.neqnV+1:end,i).^2);
        % sol_mag = sqrt(sol(1:opt.neqnV).^2+sol(opt.neqnV+1:2*opt.neqnV).^2);
        % err = norm(sol_mag-mag,inf);
        diff = mean(opt.Pr(:,i)-p_sol(xxp,yyp,study.t(i)));
        P_plot = opt.Pr(:,i)-diff;
        err = norm(sol-opt.U(:,i),1);
        % scatter3(xxp,yyp,P_plot);
        % hold on
        % scatter3(xxp,yyp,p_sol(xxp,yyp,study.t(i)));

        errp = norm(P_plot-p_sol(xxp,yyp,study.t(i)),1);
        %
        %
        if err>error(order) && i >10
            error(order) = err;%/opt.neqnV;
        end
        %
        if errp>errorp(order)  && i>10
            errorp(order) = errp;%/opt.neqnP;
        end
        % figure(101)
        % clf
        % plot(opt.U(x2_ind,i),yy(x2_ind),'-ok')
        %     [p_gridX,p_gridY] = meshgrid(xxp,yyp);
        %     contour(pGx,pGy,reshape(opt.Pr(:,100),[2*study.n_GL,2*study.n_GL]));
        %     xlim([-1,1]);ylim([-1,1]);zlim([-1,1])
        %     if mod(i,100) == 0
        %         disp('Some time')
        %     end
        %     saveas(gcf,['misc\gif\stokes_pic',num2str(i),'.png'])
        %     scatter3(xx,yy,sol1,'or')
        %
    end

    %
    % figure()
    % scatter3(xx,yy,u2,'*k')
    % hold on
    % scatter3(xx,yy,sol2,'or')
    % title('U2 velocity')
    %
    % figure()
    % scatter3(xxp,yyp,p,'*k')
    % hold on
    % % scatter3(xxp,yyp,psol,'or')
    % title('Pressure')
    % scatter3(xxp,yyp,zeros(size(psol)),'or')
    %
    % figure()
    % scatter3(xx,yy,u_mag,'*k')
    % hold on
    % scatter3(xx,yy,U_Tot,'or')
    % title('Vel Magnitude')
    % figure(101)
    % plot(opt.U(x2_ind,end),yy(x2_ind),'-','DisplayName',sprintf('Order: %d',n_GLL-1))
    %
    % [errorV] = lid_driven_benchmark(yy(x2_ind),opt.U(x2_ind,end));
    % error(order) = errorV;


    % error2(order) = norm(psol-p,inf);
    % end

end
% figure(99)
% % legend('Location','northoutside','NumColumns',4,'Interpreter','latex')
% % A = readmatrix("Roenquist_u.csv");
figure();semilogy(2*GLL-1,error,'LineWidth',3)
hold on
semilogy(2*GLL-1,errorp,'LineWidth',3)
semilogy(dbu(:,1),dbu(:,2),dbp(:,1),dbp(:,2))
% semilogy(A(:,1),A(:,2),'*r')
% grid on

% Setting the font size to 1
% set(gca, 'FontSize', 18);

% % Adding labels with LaTeX interpreter
xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);
% legend('Element 1','Element 2','Element 3','Element 4','Interpreter','latex')
% legend('Velocity','Pressure','Velocity Ron','Interpreter','latex')
legend('Velocity','Pressure','Interpreter','latex','FontSize',18)
grid on


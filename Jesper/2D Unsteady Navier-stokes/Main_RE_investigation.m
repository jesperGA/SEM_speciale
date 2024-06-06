clear
close all
clc

warning('off', 'all');

%% addpath to FEA, MESH and
addpath('SEM')
addpath('MESH')
addpath('PLOT')
addpath('misc')

% mat = [1.1,1.2,1.3,1.4;
%     1.1,1.2,1.3,1.4];
% study.p_type = 'roenquist2';
% study.p_type = 'bercover';
study.p_type = 'liddriven';
study.solve_type = 'direct'; %uzawa
if strcmp(study.solve_type,'direct') == 1

    study.direct_type = 'backslash';
    % study.direct_type = 'LU';

end
% study.solve_type = 'uzawa';
study.study_type = 'unsteady';
study.precon = 'P';
study.BC_type = 'static';
% study.BC_type = 'dynamic';

% study.study_t ype = 'steady';
if strcmp(study.study_type,'unsteady') == 1
    study.T = 0.1;


    % study.int_type = 'BDFk'; %Equivalent of solving Unsteady stokes.
    study.int_type = 'BDF3EX3'; %First order bdf for linear terms. 3 order for nonlinear terms.

    study.BDF_order = 1;

    study.U10 = 0;
    study.U20 = 0;
end

A = readmatrix("misc\lid_benchmark.txt");

%Load benchmark results
benchmark2 = A;
y_db =  benchmark2(:,1);
% px_bench = benchmark2(:,6) ;
p_db = benchmark2(:,3);
% v_bench = benchmark2(:,5);
u_db = benchmark2(:,2);
%%
n_GLL = 12;


% RE = linspace(1,2000,10);
RE = 1000;
%%

[xi,w,~] = lglnodes(n_GLL-1);
[zeta,wp] = lgwt(n_GLL-2,-1,1);
study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-2;
study.zeta = zeta;study.wp = wp;
num_el = 2;
%% MESH
[iglobV, xNV,yNV] = MeshBox_mod(1,1,num_el,num_el,n_GLL,1);
% mesh = modify_to_bercovier(xNV,yNV,iglobV);
% mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
mesh = liddriven(xNV,yNV,iglobV);
[iglobP, xNP,yNP] = MeshBox_mod(1,1,num_el,num_el,n_GLL-2,2);
% mesh_plot(mesh)

mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP,yNP];
mesh.pref_dof = 1;

plot_ind = find(xNV == 0.5);
%% [opt, study] = controller(mesh, study);

opt = [];
[opt,study] = parAssemblyQuad(mesh,opt,study);

figure(99)
hold on

% dts = logspace(-6,-1,5);
dts = 1e-6;

% p_sol = @(x,y,t) -1/4.*(cos(2.*x)+cos(2.*y)).*exp(-4.*t);


for i = 1:numel(dts)

    dt = dts(i);
    study.t = 0:dt:study.T;
    study.nt = length(study.t);
    study.dt = (study.t(2)-study.t(1));
    tic
    study.RE = RE(i);
    opt = transient_solver(opt,study,mesh);
    solve_tim = toc;
    fprintf('Time for solution @dt=%4.2e: %2.3f seconds \n',dts(i),solve_tim)


    %%

    u1 = opt.U(1:opt.neqnV,:);u2 = opt.U(opt.neqnV+1:2*opt.neqnV);p = opt.U(opt.neqnV*2+1:end);

    %% Analytical solution, validation and plots
    xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
    xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);
    fs = 18;
    for k =1:size(opt.U,2)
        figure(2)
        clf
        plot(u_db,y_db,'dr','MarkerSize',10)
        hold on
        plot(u1(plot_ind,end),yy(plot_ind),'-ok')

        % Set the x and y labels with LaTeX formatting and specified font size
        xlabel('Horizontal Velocity: $u_1(x=0.5)$', 'Interpreter', 'latex', 'FontSize', fs);
        ylabel('Position: $(x=0.5,y)$', 'Interpreter', 'latex', 'FontSize', fs);

        % Create the legend with LaTeX interpreted entries and specified font size
        legend({'Bruneau and Saad','Nodal Solution'}, 'Interpreter', 'latex', 'FontSize', fs...
            ,'Location','northoutside','NumColumns',2);

        set(gca, 'FontSize', 16);
        grid on

        % sol = opt.g_sys(xx,yy,study.t(k));
        % diff = mean(opt.Pr(:,k)-p_sol(xxp,yyp,study.t(k)));
        % P_plot = opt.Pr(:,k)-diff;
        % err = norm(sol-opt.U(:,k),1);
        % errp = norm(P_plot-p_sol(xxp,yyp,study.t(k)),1);
        % figure(101)
        % clf
        % plotSol2D(mesh,opt.U(1:opt.neqnV,i),opt.U(opt.neqnV+1:end,i),0)

    end

    % errorp(i) = errp;
    % error(i) = err;

end

% Labeling with LaTeX Interpreter
figure();loglog(dts,error,'-ok','LineWidth',3)
hold on
loglog(dts,errorp,'-oy','LineWidth',3)

xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{1}$ Error', 'Interpreter', 'latex', 'FontSize', 18);

legend('Velocity','Pressure','Interpreter','latex','FontSize',18)
grid on


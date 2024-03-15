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
study.p_type = 'roenquist2';
% study.p_type = 'bercover';
study.p_type = 'liddriven';
study.solve_type = 'direct'; %uzawa
if strcmp(study.solve_type,'direct') == 1

    % study.direct_type = 'backslash';
    study.direct_type = 'LU';

end
% study.solve_type = 'uzawa';
study.study_type = 'unsteady';
study.precon = 'P';
% study.BC_type = 'static';
study.BC_type = 'dynamic';

% study.study_t ype = 'steady';
if strcmp(study.study_type,'unsteady') == 1
    study.T = 1;


    % study.int_type = 'BDFk'; %Equivalent of solving Unsteady stokes.
    study.int_type = 'BDF1AB3'; %First order bdf for linear terms. 3 order for nonlinear terms.

    study.BDF_order = 1;

    study.U10 = 0;
    study.U20 = 0;
end

%%
n_GLL = 11;
RE = linspace(1,2000,10);
%%

[xi,w,~] = lglnodes(n_GLL-1);
[zeta,wp] = lgwt(n_GLL-2,-1,1);
study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-2;
study.zeta = zeta;study.wp = wp;
%% MESH
[iglobV, xNV,yNV] = MeshBox_mod(1,1,2,2,n_GLL,1);
% mesh = modify_to_bercovier(xNV,yNV,iglobV);
% mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
mesh = liddriven(xNV,yNV,iglobV);
[iglobP, xNP,yNP] = MeshBox_mod(1,1,2,2,n_GLL-2,2);

mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP,yNP];
mesh.pref_dof = 1;

%% [opt, study] = controller(mesh, study);

opt = [];
tic
[opt,study] = AssemblyQuad(mesh,opt,study);
dt = toc;
fprintf('Time for assembly: %2.3f seconds\n',dt)
figure(99)
hold on

dts = logspace(-6,-1,5);
dts = [1e-1 4e-2 2e-2 1e-2 5e-3 2e-3,5e-4,1e-4,5e-5];

p_sol = @(x,y,t) -1/4.*(cos(2.*x)+cos(2.*y)).*exp(-4.*t);


for i = 1:numel(dts)

    dt = dts(i);
    study.t = 0:dt:study.T;
    study.nt = length(study.t);
    study.dt = (study.t(2)-study.t(1));
    tic
    study.RE = 1;
    opt = transient_solver(opt,study,mesh);
    solve_tim = toc;
    fprintf('Time for solution @dt=%4.2e: %2.3f seconds \n',dts(i),solve_tim)


    %%

    u1 = opt.U(1:opt.neqnV);u2 = opt.U(opt.neqnV+1:2*opt.neqnV);p = opt.U(opt.neqnV*2+1:end);

    %% Analytical solution, validation and plots
    xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
    xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);

    for k =1:size(opt.U,2)
        sol = opt.g_sys(xx,yy,study.t(k));
        diff = mean(opt.Pr(:,k)-p_sol(xxp,yyp,study.t(k)));
        P_plot = opt.Pr(:,k)-diff;
        err = norm(sol-opt.U(:,k),1);
        errp = norm(P_plot-p_sol(xxp,yyp,study.t(k)),1);
    end

    errorp(i) = errp;
    error(i) = err;

end

% Labeling with LaTeX Interpreter
figure();loglog(dts,error,'-ok','LineWidth',3)
hold on
loglog(dts,errorp,'-oy','LineWidth',3)

xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{1}$ Error', 'Interpreter', 'latex', 'FontSize', 18);

legend('Velocity','Pressure','Interpreter','latex','FontSize',18)
grid on


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
% study.p_type = 'roenquist';
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
% study.study_t ype = 'steady';
if strcmp(study.study_type,'unsteady') == 1
    study.T = 0.2;
    study.nt = 2000;
    study.t = linspace(0,study.T,study.nt);
    study.dt = (study.t(2)-study.t(1));

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

for i = 1:numel(RE)
    tic
    study.RE = RE(i);
    opt = transient_solver(opt,study,mesh);
    dt = toc;
    fprintf('Time for solution @RE=%4.2f: %2.3f seconds',RE(i),dt)


    %%

    u1 = opt.U(1:opt.neqnV);u2 = opt.U(opt.neqnV+1:2*opt.neqnV);p = opt.U(opt.neqnV*2+1:end);

    %% Analytical solution, validation and plots
    xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
    x2_ind = find(xx==0.5);
    xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);


    figure(99)
    plot(opt.U(x2_ind,end),yy(x2_ind),'-o');

end
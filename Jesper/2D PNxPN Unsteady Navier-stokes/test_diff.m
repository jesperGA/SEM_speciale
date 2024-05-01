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
% study.p_type = 'liddriven';
study.solve_type = 'direct'; %uzawa
if strcmp(study.solve_type,'direct') == 1

    % study.direct_type = 'backslash';
    study.direct_type = 'LU';

end
study.P_order = 'PnPn-2';
study.P_order = 'PnPn';

% study.solve_type = 'uzawa';
study.study_type = 'unsteady';
% study.precon = 'mhat';
study.precon = 'P';
% study.study_t ype = 'steady';
if strcmp(study.study_type,'unsteady') == 1
    study.T = 1;
    % study.nt = 1e3;
    % study.t = linspace(0,study.T,study.nt);
    study.t = 0:1e-4:study.T;
    study.nt = length(study.t);
    study.dt = (study.t(2)-study.t(1));

    % study.int_type = 'BDFk'; %Equivalent of solving Unsteady stokes.
    study.int_type = 'BDF1AB3'; %First order bdf for linear terms. 3 order for nonlinear terms.
    study.RE = 1;
    study.BDF_order = 1;

    study.U10 = 0;
    study.U20 = 0;

    study.BC_type = 'dynamic';
    % study.BC_type = 'static';
end
% % GLL = 5:1:14;
n_GLL = 4;
order = 1;

error(order) = 0;
errorp(order) = 0;

[xi,w,~] = lglnodes(n_GLL-1);
[zeta,wp] = lgwt(n_GLL-2,-1,1);
study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-2;
study.zeta = zeta;study.wp = wp;
%% MESH
[iglobV, xNV,yNV] = MeshBox_mod(2,2,2,2,n_GLL,1);
% mesh = modify_to_bercovier(xNV,yNV,iglobV);
% mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
mesh = modify_to_roenquist_mesh(xNV,yNV,iglobV);
% mesh = liddriven(xNV,yNV,iglobV);
% mesh = displaceInternalNodes(mesh,0,-0.5);
[iglobP, xNP,yNP] = MeshBox_mod(2,2,2,2,study.n_GL,2);

mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP-1,yNP-1];
mesh.pref_dof = 1;
% mesh_plot(mesh);
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
% figure(101)
% db = readmatrix("misc\lid_benchmark.txt");
% plot(db(:,2),db(:,1),'dk','DisplayName','Benchmark')
% hold on

opt = [];
[opt,study] = PnPn_AssemblyQuad(mesh,opt,study);
%%%%%%%%%%%%%%%%%%%%% SLET
U = opt.g_sys(mesh.Xv(:,2),mesh.Xv(:,3),0);
u2 = U(opt.neqnV+1:end);
u1 = U(1:opt.neqnV);

du1dx1_fisse = opt.DE1 * u1;
du1dx2_fisse = opt.DE2 * u1;
du2dx1_fisse = opt.DE1 * u2;
du2dx2_fisse = opt.DE2 * u2;

du1dx1 = @(X1,X2)  sin(X1).*sin(X2).*exp(-2*0);
du1dx2 = @(X1,X2) -cos(X1).*cos(X2).*exp(-2*0);
du2dx1 = @(X1,X2)  cos(X1).*cos(X2).*exp(-2*0);
du2dx2 = @(X1,X2) -sin(X1).*sin(X2).*exp(-2*0);

du1dx1_pik = du1dx1(mesh.Xv(:,2),mesh.Xv(:,3));
du1dx2_pik = du1dx2(mesh.Xv(:,2),mesh.Xv(:,3));
du2dx1_pik = du2dx1(mesh.Xv(:,2),mesh.Xv(:,3));
du2dx2_pik = du2dx2(mesh.Xv(:,2),mesh.Xv(:,3));

figure()
scatter3(mesh.Xv(:,2),mesh.Xv(:,3),du2dx2_fisse)
hold on
scatter3(mesh.Xv(:,2),mesh.Xv(:,3),du2dx2_pik)
%%%%%%%%%%%%%%%%%%%%%%% SLET SLUT

%%%%%%%%%%%%%%%%%%%%% SLET
% u1 = study.U1(mesh.Xv(:,2),mesh.Xv(:,3),0);
% u2 = study.U2(mesh.Xv(:,2),mesh.Xv(:,3),0);

ddu1ddx2_fisse = opt.DD1* u2;
ddu2dx1dx2_fisse = opt.DE1 * (opt.DE2 * u2);
ddu2ddx1_fisse = opt.DD2*u2;
ddu1dx1dx2_fisse = opt.DE1 * (opt.DE2 * u1);

ddu1ddx2   = @(X1,X2)  cos(X1).*sin(X2).*exp(-2*0);
ddu2dx1dx2 = @(X1,X2) -cos(X1).*sin(X2).*exp(-2*0);
ddu2ddx1   = @(X1,X2) -sin(X1).*cos(X2).*exp(-2*0);
ddu1dx1dx2 = @(X1,X2)  sin(X1).*cos(X2).*exp(-2*0);

ddu1ddx2_pik   = ddu1ddx2(mesh.Xv(:,2),mesh.Xv(:,3));
ddu2dx1dx2_pik = ddu2dx1dx2(mesh.Xv(:,2),mesh.Xv(:,3));
ddu2ddx1_pik   = ddu2ddx1(mesh.Xv(:,2),mesh.Xv(:,3));
ddu1dx1dx2_pik = ddu1dx1dx2(mesh.Xv(:,2),mesh.Xv(:,3));

figure()
scatter3(mesh.Xv(:,2),mesh.Xv(:,3),ddu2ddx1_fisse)
hold on
scatter3(mesh.Xv(:,2),mesh.Xv(:,3),ddu2ddx1_pik)
%%%%%%%%%%%%%%%%%%%%%%% SLET SLUT
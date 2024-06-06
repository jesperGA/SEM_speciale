% clear
close all
clc

warning('off', 'all');

%% addpath to FEA, MESH and
addpath('SEM')
addpath('MESH')
addpath('PLOT')
addpath('misc')

study = study_settings;
% mat = [1.1,1.2,1.3,1.4;
%     1.1,1.2,1.3,1.4];

% GLL = 3:1:10;
GLL = 2;
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
nelx = 32; nely = 16;
h = 1; L = 2;
V0 = 1;
element_to_remove = 5;
N_interp = 5;

for order = 1:numel(GLL)
    n_GLL = GLL(order);

    error(order) = 0;
    errorp(order) = 0;

    [xi,w,~] = lglnodes(n_GLL-1);
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-1;
    [zeta,wp] = lgwt(study.n_GL,-1,1);
    study.zeta = zeta;study.wp = wp;
    %% MESH
    [iglobV, xNV,yNV] = MeshBox_mod(L,h,nelx,nely,n_GLL,1);
    [iglobP, xNP,yNP] = MeshBox_mod(L,h,nelx,nely,study.n_GL,2);
    % [iglobINT, xNINT,yNINT] = MeshBox_mod(L,h,nelx,nely,5,3);

    mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP,yNP];
    mesh.IXv = iglobV;mesh.Xv = [(1:numel(xNV)).',xNV,yNV];

    if isfield(study,'Brinkmann') && study.Brinkmann
        mesh = modify_to_pipe_brink(h,L,mesh,V0,1);
    else
        mesh = modify_to_Backstep(h,L,mesh,V0,element_to_remove,1);
    end
    % mesh.pref_dof = 1;

    %% SOLVE SYSTEM
    [opt, study] = controller(mesh, study);

    %% POST PROC
    % X_int = [(1:numel(xNINT))',xNINT,yNINT];
    % [X_int,iglobINT,~] = remove_element(X_int,iglobINT,element_to_remove);
    % mesh_int.IX = iglobINT; mesh_int.X = X_int;
    name = sprintf('backstep_Order%dRE%d',n_GLL-1,study.RE);
    save2nek5000(name,opt.movie_step,opt,mesh,study,N_interp)

    u1 = opt.U(1:opt.neqnV);u2 = opt.U(opt.neqnV+1:2*opt.neqnV);p = opt.U(opt.neqnV*2+1:end);

    %% Analytical solution, validation and plots
    xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
    xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);

    % Initialize the GIF file
    filename = 'output_animation.gif';  % Define the filename for the GIF

    % for i = 1:size(opt.U,2)
    %     figure(100)
    %     clf
    %     % Uncomment the following lines if you want to use scatter plots
    %     % scatter3(xx, yy, opt.U(1:opt.neqnV,i), '*k');
    %     % scatter3(xxp, yyp, opt.Pr(:,i), '*k');
    % 
    %     % Plotting solution - make sure you replace 'mesh' and other variables with your actual data
    %     plotSol2D(mesh, opt.U(1:opt.neqnV,i), opt.U(opt.neqnV+1:end,i), 0);
    % 
    %     % Capture the frame
    %     % drawnow;
    %     % frame = getframe(100);
    %     % im = frame2im(frame);
    %     % [imind, cm] = rgb2ind(im, 256);
    %     % 
    %     % Write to the GIF File
    %     % if i == 1
    %     %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.03);
    %     % else
    %     %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
    %     % end
    % end



end
% figure(99)
% % legend('Location','northoutside','NumColumns',4,'Interpreter','latex')
% % A = readmatrix("Roenquist_u.csv");
% figure();semilogy(2*GLL-1,error,'LineWidth',3)
% hold on
% semilogy(2*GLL-1,errorp,'LineWidth',3)
% semilogy(dbu(:,1),dbu(:,2),dbp(:,1),dbp(:,2))
% semilogy(A(:,1),A(:,2),'*r')
% grid on

% Setting the font size to 1
% set(gca, 'FontSize', 18);

% % Adding labels with LaTeX interpreter
% xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);
% legend('Element 1','Element 2','Element 3','Element 4','Interpreter','latex')
% legend('Velocity','Pressure','Velocity Ron','Interpreter','latex')
% legend('Velocity','Pressure','Interpreter','latex','FontSize',18)
% grid on


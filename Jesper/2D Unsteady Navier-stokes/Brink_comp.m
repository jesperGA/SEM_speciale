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
GLL = 14;
% n_interp = 20;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
nelx = 10; nely = 2;
h = 1; L = 5;
V0 = 1;
element_to_remove = 5;
N_interp = 20;

n_GLL = GLL(order);

error(order) = 0;
errorp(order) = 0;

[xi,w,~] = lglnodes(n_GLL-1);
[zeta,wp] = lgwt(n_GLL-2,-1,1);
study.xi = xi;study.w = w;study.n_GLL = n_GLL;study.n_GL = n_GLL-2;
study.zeta = zeta;study.wp = wp;
%% MESH
[iglobV, xNV,yNV] = MeshBox_mod(L,h,nelx,nely,n_GLL,1);
[iglobP, xNP,yNP] = MeshBox_mod(L,h,nelx,nely,study.n_GL,2);
% [iglobINT, xNINT,yNINT] = MeshBox_mod(L,h,nelx,nely,5,3);

mesh.IXp = iglobP;mesh.Xp = [(1:numel(xNP)).',xNP,yNP];
mesh.IXv = iglobV;mesh.Xv = [(1:numel(xNV)).',xNV,yNV];

for i = 1:2
    if i==2
        study.Brinkmann = 0;
    end
    if isfield(study,'Brinkmann') && study.Brinkmann
        mesh = modify_to_pipe_brink(h,L,mesh,V0);
    else
        mesh = modify_to_pipe(h,L,mesh,V0,element_to_remove);
    end
    % mesh.pref_dof = 1;

    %% SOLVE SYSTEM
    [opt, study] = controller(mesh, study);

    %% POST PROC
    % X_int = [(1:numel(xNINT))',xNINT,yNINT];
    % [X_int,iglobINT,~] = remove_element(X_int,iglobINT,element_to_remove);
    % mesh_int.IX = iglobINT; mesh_int.X = X_int;
    % name = sprintf('pipe_whole_brinkmanRE%d',study.RE);
    % save2nek5000(name,1:200:study.nt,opt,mesh,study,N_interp)
    [A] = twoD_element_interpolator(mesh,51,opt.U(1:opt.neqnV,end),opt.U(opt.neqnV+1:end,end),opt.Pr(:,end));
    ind_mat = A(:,:,2)==0.3;
    
    u1 = opt.U(1:opt.neqnV,:);u2 = opt.U(opt.neqnV+1:2*opt.neqnV,:);p = opt.Pr;

    %% Interpolated values
    xx = A(:,:,1);
    xx = xx(ind_mat) ;
    [xx,sort_ind] = sort(xx);

    yy = A(:,:,2);
    yy = yy(ind_mat);

    uu = A(:,:,3);
    uu = uu(ind_mat);

    vv = A(:,:,4);
    vv = vv(ind_mat);

    pp = A(:,:,5);
    pp = pp(ind_mat);
    pp = pp(sort_ind);

    vel_mag = uu.^2+vv.^2;
    vel_mag = vel_mag(sort_ind);

    figure(1);
    plot(xx,vel_mag,'-o')
    hold on

    figure(2);
    plot(xx,pp,'-o');
    hold on

    data{i} = [xx,vel_mag,pp];

end

for i =1:2

    figure(i)
    xlabel('$x$-position', 'Interpreter', 'latex', 'FontSize', 20);
    if i ==1
    ylabel('Velocity @ $y=0.3$', 'Interpreter', 'latex', 'FontSize', 20);
    else
        ylabel('Pressure @ $y=0.3$', 'Interpreter', 'latex', 'FontSize', 20);
    end
    legend('Brinkmann term','Element removal','Interpreter','latex','FontSize',20)
    set(gca, 'FontSize', 18);
end
% Initialize the GIF file
% filename = 'output_animation.gif';  % Define the filename for the GIF

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


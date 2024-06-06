addpath('SEM')
addpath('MESH')
addpath('PLOT')
addpath('misc')

close all
clear all

files = dir('N*.mat');
files = files(~cellfun(@isempty, regexp({files.name}, '^N\d+\.mat$')));
N_interp = 50;
pos = 0.5;

subfolder = 'pics_order_comp';
% Check if the subfolder exists, if not, create it
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

for i = 1:4
    file_name = files(i).name;
    load(file_name);

    N = study.n_GLL-1;

    ind = mesh.Xv(:,3) == 0.5;
    x = mesh.Xv(ind,2);
    u = opt.U(ind,end);

    tempv = opt.U(opt.neqnV+1:end,end);
    v = tempv(ind);
    u_mag2 = sqrt(u.^2+v.^2);
    % u_mag2 = v;


    [A] = twoD_element_interpolator(mesh,N_interp,opt.U(1:opt.neqnV,end),opt.U((opt.neqnV+1):end,end),opt.Pr(:,end));

    ind_mat = A(:,:,2)==0.5;

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
    
    % vel_mag = vv;
    vel_mag = sqrt(uu.^2+vv.^2);
    vel_mag = vel_mag(sort_ind);

    % figure(i);
    figure('position',[500,500,900,450])
    plot(xx,vel_mag,'-k','LineWidth',2)
    hold on
    plot(x,u_mag2,'kd')
    xlabel('$x$-position', 'Interpreter', 'latex', 'FontSize', 24);

    ylabel('Velocity @ $y=0.5$', 'Interpreter', 'latex', 'FontSize', 24);
    legend('Interpolated values','Nodal value','Interpreter','latex','FontSize',24,'Location','northoutside','Orientation','horizontal')
    set(gca, 'FontSize', 20);
    axis equal
    ylim([0,2.8]);xlim([0,5])

    filename = sprintf('v_order_comp%d.pdf', N);  % Using PNG format, you can change the format as needed
    % Full file path
    fullFilePath = fullfile(subfolder, filename);
    % Save the current figure to the specified path
    exportgraphics(gca, fullFilePath, 'ContentType', 'vector', 'BackgroundColor', 'none');
    % exportgraphics(gca, fullFilePath, 'BackgroundColor', 'none');

    % Make sure the figure is the current figure
    % figure(gcf);

    % Save the figure
    % print(gcf, fullFilePath, '-dsvg', '-r300');

    if N==15
        ylim([1.4,2.1]);xlim([1.4,2])
        filename = sprintf('v_order_comp%d_zoom.svg', N);  % Using PNG format, you can change the format as needed
        % Full file path
        fullFilePath = fullfile(subfolder, filename);
        % Save the current figure to the specified path
        % exportgraphics(gca, fullFilePath, 'BackgroundColor', 'none');

        % Make sure the figure is the current figure
        figure(gcf);

        % Save the figure
        % print(gcf, fullFilePath, '-dsvg', '-r300');
    end


    % figure(2);
    % plot(xx,pp,'-k','LineWidth',2);
    % hold on

end



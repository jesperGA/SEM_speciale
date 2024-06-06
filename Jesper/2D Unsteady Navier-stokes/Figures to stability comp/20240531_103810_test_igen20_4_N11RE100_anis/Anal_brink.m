load("allData.mat")
subfolder = 'pics_order_comp';
% Check if the subfolder exists, if not, create it
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

N_interp = 25;
pos = 0.205;

N = study.n_GLL-1;

ind = mesh.Xv(:,3) == pos;
x = mesh.Xv(ind,2);
u = opt.U(ind,end);

tempv = opt.U(opt.neqnV+1:end,end);
v = tempv(ind);
u_mag2 = sqrt(u.^2+v.^2);
% u_mag2 = v;


[A] = twoD_element_interpolator(mesh,N_interp,opt.U(1:opt.neqnV,end),opt.U((opt.neqnV+1):end,end),opt.Pr(:,end));

ind_mat = A(:,:,2)==pos;

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

ylabel('Velocity @ $y=h/2$', 'Interpreter', 'latex', 'FontSize', 24);
legend('Interpolated values','Nodal value','Interpreter','latex','FontSize',24,'Location','northoutside','Orientation','horizontal')
set(gca, 'FontSize', 20);
% axis equal
ylim([0,2]);xlim([0,2.05])
fill([0.15,0.25,0.25,0.15],[1e-16,1e-16,1e3,1e3],[0.8,0.8,0.8],'FaceAlpha',0.3,...
    'HandleVisibility','off');

filename = sprintf('v_order_comp%d.pdf', N);  % Using PNG format, you can change the format as needed
% Full file path
fullFilePath = fullfile(subfolder, filename);
% Save the current figure to the specified path
exportgraphics(gca, fullFilePath, 'ContentType', 'vector', 'BackgroundColor', 'none');
% exportgraphics(gca, fullFilePath, 'BackgroundColor', 'none');

figure('position',[500,500,900,450])
yyaxis left
semilogy(xx,vel_mag,'-k','LineWidth',2)
hold on
plot(x,u_mag2,'kd')
fill([0.15,0.25,0.25,0.15],[1e-16,1e-16,1e3,1e3],[0.8,0.8,0.8],'FaceAlpha',0.3,...
    'HandleVisibility','off');
ylim([1e-3,0.2])
ylabel('Velocity @ $y=h/2$', 'Interpreter', 'latex', 'FontSize', 24);
% set(gca, 'YScale', 'log')
set(gca, 'YColor', 'black'); 
set(gca, 'FontSize', 20);

yyaxis right

temp = full(diag(opt.alpha));
plot(x,temp(ind),'-ob','LineWidth',2)
ylim([0,0.15])
xlim([0.1,0.3])
set(gca, 'YColor', 'blue'); 
set(gca, 'FontSize', 20);

xlabel('$x$-position', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Brinkmann term: $\alpha$', 'Interpreter', 'latex', 'FontSize', 24);

filename = 'zoom.pdf';
fullFilePath = fullfile(subfolder, filename);
exportgraphics(gca, fullFilePath, 'ContentType', 'vector', 'BackgroundColor', 'none');

% Find linje i midten af mesh
% og plot filter med R som ticks, mpsje bare brink langs x også....
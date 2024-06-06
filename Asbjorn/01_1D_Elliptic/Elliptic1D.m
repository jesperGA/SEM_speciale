clear 
close all
clc

% Add necessary paths
addpath('FEM')
addpath('MESH')
addpath('PLOT')
addpath('TOPOPT')
addpath('Saved_values')

study.maxit=1; 
study.p = 1;

% Call the mesh routine 
L = pi;
ne = 2;
for i =1:15
    study.N = i;
    if study.N==1
        study.etype = 'Linear';
    elseif study.N>1
        study.etype = 'SEM';
    end
    
    [mesh] = mesh1D(L, ne, study.N);
    
    % % Plot the mesh
    % fig=figure();
    % plot_mesh(mesh);
    % enhance_plot(0,0,0,0,0)
    % fig.Position = [744 958.6000 560 119.4000];
    % % saveas(gcf,'mu0','epsc')
    
    study.f = exp(mesh.X(:,2)) .* (cos(mesh.X(:,2))-sin(mesh.X(:,2)));
    % study.f = (cos(mesh.X(:,2))-sin(mesh.X(:,2)));
    
    % % Plot f
    % defaultColors = get(groot,'DefaultAxesColorOrder');
    % fig=figure();
    % subplot(2,1,1)
    % plot(mesh.X(:,2), study.f,'Color', defaultColors(1, :))
    % xlabel('x')
    % ylabel('f(x)')
    % % xlims2 = [time*(-2/10)/T_it,time*(2/10)/T_it];
    % % xlim([xlims2])
    % enhance_plot(0, 0, 0, 0, 0)
    
    opt = Controller(mesh, study);
    
    % plot_timeseries(mesh,opt_empty,study.n_t)
    x=mesh.X(:,2);
    solution = -exp(x).*cos(x)/2 - exp(x).*sin(x)/2 - ((1 + exp(pi)).*x)./(2*pi) + 1/2;
    solution = -sin(x);
    e=abs(opt.U-solution);
    norme(i) = norm(e,'inf');
    dofs(i) = study.N*2+1;
end

% figure()
% semilogy(dofs,norme,'o')

% Mindste kvadraters metode
x=[dofs(3:end)]';
y=[log10(norme(3:end))]';
A=[x.^0,dofs(3:end)'];
c = (A'*A)\(A'*y);


% hold on
% semilogy(dofs(3:end),10.^(c(2).*dofs(3:end)+c(1)))
% enhance_plot(0,0,0,0,0)



F = @(x) (-sin(x));
plotNodalSolution(F,mesh.X,opt.U)

p1 = polyfit((dofs(3:end)),log10(norme(3:end)),1);

% 
% figure()
% semilogy(dofs,norme,'o')
% hold off
% grid on

% for k=5:11
% q=log((norme(k+1)-norme(k))/(norme(k)-norme(k-1)))/log((norme(k)-norme(k-1))/(norme(k-1)-norme(k-2)));
% end
% 
% for i=3:11
%     seq(i)=norme(i+1)/norme(i);
% end
% plot(seq)
% ylim([0 1])

plotConvergence(norme,dofs)

plotMeshh(mesh)

% figure
% hold on
% for i=1:2
% plot(opt.elementX(i,:),opt.elementSolution(i,:))
% end

% Plot the sparsity pattern of the matrix
figure;
spy(opt.A_noBCs, 25); % 'k' specifies the color black, and 25 is the marker size
% title('Sparsity Pattern of Matrix A');
% xlabel('DoF Index');
% ylabel('DoF Index');
axis off
enhance_plot(0, 22, 0, 15, 0);
saveas(gcf,'Sparcity_pattern_A','epsc')

% Plot the sparsity pattern of the matrix
figure;
spy(opt.B_noBCs, 25); % 'k' specifies the color black, and 25 is the marker size
% title('Sparsity Pattern of Matrix B');
% xlabel('DoF Index');
% ylabel('DoF Index');
axis off
enhance_plot(0, 22, 0, 15, 0);
saveas(gcf,'Sparcity_pattern_B','epsc')

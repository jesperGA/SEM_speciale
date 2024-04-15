clear
close all
clc
mat = [1.1,1.2,1.3,1.4;
    1.1,1.2,1.3,1.4];

GLL = 3:1:20;
% GLL = 1;
for i = 1:numel(GLL)
    % i = 1;
    n_GLL = GLL(i); %Specify number of GLL points
    % n_GLL = 6;

    L = pi;

    [xi,w,~] = lglnodes(n_GLL-1);
    % xi = 
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;
    mesh = [];
    mesh = regular_bragg_grating(L,3,n_GLL,mat,study);
    % mesh.X(:,2) = linspace(mesh.X(1,2),mesh.X(end,2),length(mesh.X));
    
    % study.xi = linspace(-1,1,n_GLL);
    % study.w(:) = sum(w)/length(study.xi);
    % figure(1)
    % plot(mesh.X(:,2),ones(size(mesh.X,1)),'ok');

    opt = [];

    [opt,study] = Assembly1bar(mesh,opt,study);
    xv = mesh.X(:,2);
    f = exp(xv).*(cos(xv)-sin(xv));
    A = opt.K;B = opt.M;
    f(1) = 0;f(end) = 0;
    A(1,:) = 0;A(:,1) = 0;
    A(end,:) = 0;A(:,end) = 0;
    A(1) = 1; A(end) = 1;


    u = A\(B*f);

    xx = linspace(0,pi,200);
    analytical = -sin(xx);

    sol = -sin(xv);
    
    [data] = oneD_element_interpolator(mesh,u);
    analytical2 = -sin(data(:,1));
    if i == 1|| i==4 || i==15 || i==18
        figure()
        % analytical = -exp(xv).*cos(xv)/2 - exp(xv).*sin(xv)/2 - ((1 + exp(pi)).*xv)./(2*pi) + 1/2;


        plot(xv,u,'ok','LineWidth',2);
        hold on

        plot(xx,analytical,data(:,1),data(:,2),'LineWidth',2)


        legend("Numerical","analytical","Interpolated data")
        title(['order = ',num2str(n_GLL-1)])
    end

    err(i) = norm(sol - u,2);
    err2(i) = norm(analytical2-data(:,2),2);



end

figure();semilogy(2*GLL+1,err,'-ok','LineWidth',3)
hold on
semilogy(2*GLL+1,err2,'-or','LineWidth',3)
grid on

% Setting the font size to 18
set(gca, 'FontSize', 18);

% Adding labels with LaTeX interpreter
xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);

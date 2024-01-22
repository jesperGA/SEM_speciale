clear
close all
clc
mat = [1.1,1.2,1.3,1.4;
    1e6,1.2,1.3,1.4];
% mat = ones([2,10]);
% mat(2,:) = linspace(2,1e9,10);
scale = linspace(2,1e9,10);
GLL = 3:10:80;
% GLL = 1;
for i = 1:length(GLL)
    % i = 1;
    n_GLL = GLL(i); %Specify number of GLL points
    % n_GLL = 10;

    L = pi;

    [xi,w,~] = lglnodes(n_GLL-1);
    % xi =
    study.xi = xi;study.w = w;study.n_GLL = n_GLL;
    mesh = [];
    mesh = regular_bragg_grating(L,3,n_GLL,mat,study);
    mesh.material = mat(:,1);
    % figure(1)
    % plot(mesh.X(:,2),ones(size(mesh.X,1)),'ok');

    opt = [];

    [opt,study] = Assembly1bar(mesh,opt,study);
    xv = mesh.X(:,2);
    f = ones(length(xv),1)*pi;
    quart = floor(1/4*length(xv));
    f(quart:end-quart) = f(quart:end-quart)*scale(end);




    A = opt.K;B = opt.M;
    f(1) = 0;f(end) = 0;
    A(1,:) = 0;A(:,1) = 0;
    A(end,:) = 0;A(:,end) = 0;
    A(1) = 1; A(end) = 1;


    u = A\(B*f);

    xx = linspace(0,pi,200);
    analytical = -sin(xx);

    [data] = oneD_element_interpolator(mesh,u);
    analytical2 = -sin(data(:,1));

    figure()
    % analytical = -exp(xv).*cos(xv)/2 - exp(xv).*sin(xv)/2 - ((1 + exp(pi)).*xv)./(2*pi) + 1/2;


    plot(xv,u,'ok','LineWidth',2);
    hold on

    plot(data(:,1),data(:,2),'LineWidth',2)


    legend("Numerical","Interpolated data")
    title(['N = ',num2str(GLL(i))])




end

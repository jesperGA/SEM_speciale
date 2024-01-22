clear 
close all
clc

% Add necessary paths
addpath('FEM')
addpath('MESH')
addpath('PLOT')
addpath('TOPOPT')
addpath('Saved_values')

N=5;
[xi,w]=lglnodes(N);

% Cardinal
figure()
for e = 1:N+1
    plot(xi,zeros(N+1,1),'o')
    hold on
    plot(linspace(-1,1),CardinalPolynomial(xi,e,linspace(-1,1)),'Color','red')
    enhance_plot(0,0,0,0,0)
    hold off
end


% Legendre
xi=linspace(-1,1);

figure()
plot(lglnodes(N),zeros(N+1,1),'o')
hold on
xticks(lglnodes(N))
grid on
for i=1:100
    L_N(i) = 1;
end
plot(xi,L_N)
enhance_plot(0,0,0,0,0)
hold off
for N=1:5
    plot(lglnodes(N),zeros(N+1,1),'o')
    hold on
    xticks(lglnodes(N))
    grid on
    for i=1:100
        L_N(i) = LegendrePoly(N,xi(i));
    end
    plot(xi,L_N)
    enhance_plot(0,0,0,0,0)
    hold off
end

function L_N = LegendrePoly(N,xi)
    L=zeros(N+1,1);
    L(1) = 1;
    L(2) = xi;
    for i=3:N+1
        n=i-1;
        L(i) = 1/n*((2*n-1)*xi*L(i-1)-(n-1)*L(i-2));
    end
    L_N=L(end);
end
clear 
close all
clc

% Add necessary paths
addpath('FEM')
addpath('MESH')
addpath('PLOT')

%---------------------

clear all; close all;clc
x1lim=[-1 1];
x2lim=[-1 1];
F= @(X,Y) 1/4.*(1-X).*(1-Y)
for i=1:length(X)
    for j= 1:length(Y)
        lk = CardinalPolynomial(lglnodes(99),k,xi(j))
    end
end


figure(1),
x1=linspace(x1lim(1),x1lim(2),1e2);
x2=linspace(x2lim(1),x2lim(2),1e2);
[X,Y]=meshgrid(x1,x2);
mesh(X,Y,F(X,Y)), grid minor,
% xlabel('x_1','fontsize',14),
% ylabel('x_2','fontsize',14),
% zlabel('f (x)','fontsize',14),
% hold on
% plot3(1,1,0,'r*')
% enhance_plot(0,0,0,0,0)
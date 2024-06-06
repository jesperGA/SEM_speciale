% Clean up the workspace and figures to start fresh
clear; close all; clc;

% Add necessary paths for project dependencies
addpath('FEM', 'MESH', 'PLOT');

defaultColors = get(groot, 'DefaultAxesColorOrder');

% Initialize variables
N = 7; % Degree of polynomial
[xi, w] = lglnodes(N); % Legendre-Gauss-Lobatto nodes and weights
[xi_gl,w_gl]=lgwt(N-2+1,-1,1); 

    % % Compute the corresponding function values
    % y = f(xi);
[xi,w,h]=GetGLL(N+1);
z=linspace(-1,1,1e2);
% z=xi
plotto = 2;

fntsize = 22;
% z = linspace(-1, 1, 100); % Define z if not already defined

% Calculate Legendre Polynomials
L = zeros(length(z), 10); % Initialize L with zeros
for NN = 0:9
    for i = 1:length(z)
        L(i,NN+1) = legendreP(NN, z(i)); % Use legendreP function
    end
end

% Plotting Even Legendre Polynomials L0, L2, ..., L8
figure();
plot(z, L(:,1:2:end), 'LineWidth', 1.5); % Adjust as needed
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize);
ylabel('$L$', 'Interpreter', 'latex', 'FontSize', fntsize);
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
grid on;
legend(arrayfun(@(n) ['$L_{' num2str(n) '}$'], 0:2:8, 'UniformOutput', false), ...
       'Location', 'eastoutside', 'Interpreter', 'latex'); % Legends in a column to the right
enhance_plot(0, fntsize, 0, 0, 0);
% saveas(gcf,'Figures/L_even','epsc');

% Plotting Odd Legendre Polynomials L1, L3, ..., L9
figure();
plot(z, L(:,2:2:end), 'LineWidth', 1.5); % Adjust as needed
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize);
ylabel('$L$', 'Interpreter', 'latex', 'FontSize', fntsize);
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
grid on;
legend(arrayfun(@(n) ['$L_{' num2str(n) '}$'], 1:2:9, 'UniformOutput', false), ...
       'Location', 'eastoutside', 'Interpreter', 'latex'); % Legends in a column to the right
enhance_plot(0, fntsize, 0, 0, 0);
% saveas(gcf,'Figures/L_odd','epsc');

n = N + 1; % Increment N by 1
fig = figure(); % Create a new figure
fig.Position = [744.0000  750.6000  667.4000  327.4000];
hold on
zz=(n ) .* (L(:, n-1) - z' .* L(:, n));
plot(z, zz); % Plot the main function
xlim([-1.1, 1.1]); % Set the x-axis limits
ylim([min(zz)-1 max(zz)+1])
xticks(xi); % Set the x-ticks
grid on; % Enable the grid
% ylim([-1.1, 1.1]); % Set the y-axis limits if needed
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
ylabel(['$\left(1-\xi^2\right) L_', num2str(N), '^{\prime}(\xi)$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting

yline(0, 'Color', 'k', 'LineStyle', '-')
plot(xi,zeros(length(xi),1),'o')
legend('','','GLL Points (Roots)', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside')
enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
xtickformat('%.2f'); % Set x-tick labels format
% saveas(gcf,['Figures\GLL_Points_N',num2str(N)],'epsc')

plot_func = 3;
for N =[3 7]
    for j=[2 1]
        distributions={'GLL','lin'};
        dist = distributions{j};
        if strcmp(dist,'GLL')
            xi = lglnodes(N);
        elseif strcmp(dist,'lin')
            xi = linspace(-1,1,N+1)';
        end
        for j = 1:N+1
            for i = 1:length(z)
                pi(i,j) = Lagrange(j, xi', z(i));
                D(i,j) = derivativeLagrange(j, xi', z(i));
            end
        end
        
        % Plot Lagrange
        figure;
        plot(z, pi(:,plot_func)); % Plot cardinal polynomial
        hold on;
        plot([xi; xi(plot_func)], [zeros(N+1, 1); 1], 'o'); % Plot nodes
        xlim([-1.1, 1.1]); % Set the x-axis limits
        % if N==3
        %     ylim([-0.35 1.1])
        % elseif N==7
        %     ylim([-1.65 1.15])
        % end
        ylim([-0.4 1.1])
        % xticks(xi); % Set the x-ticks
        grid on; % Enable the grid
        % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
        xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
        ylabel(['$\left(1-\xi^2\right) L_', num2str(N), '^{\prime}(\xi)$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
        % legend('','','GLL Points (Roots)', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside')
        axis off
        enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
        % xtickformat('%.2f'); % Set x-tick labels format
        line([xi(plot_func) xi(plot_func)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line
        line([[xi(1) xi(end)]],[0 0], 'Color', defaultColors(1,:), 'LineStyle', '-','LineWidth',0.7)
        saveas(gcf,['Figures\point_dist_N_',num2str(N),dist],'epsc')

    end
end

for N=[1:3 5]
    % if N==1
    %     set(groot, 'defaultAxesXColor', defaultColors(1,:), 'defaultAxesYColor', defaultColors(1,:), 'defaultAxesZColor', defaultColors(1,:));
    % elseif N==2
    %     set(groot, 'defaultAxesXColor', defaultColors(2,:), 'defaultAxesYColor', defaultColors(2,:), 'defaultAxesZColor', defaultColors(2,:));
    % elseif N>2
    %     set(groot, 'defaultAxesXColor', defaultColors(4,:), 'defaultAxesYColor', defaultColors(4,:), 'defaultAxesZColor', defaultColors(4,:));
    % end
fig=figure;
fig.Position=[746.6000 526.6000 341.6000 245.6000]
hold on
xi = lglnodes(N);
plot([xi], [zeros(N+1, 1)], 'o', 'Color', defaultColors(1,:),'Color','k'); % Plot nodes
plot([xi], [ones(N+1, 1)], 'o', 'Color', defaultColors(1,:),'Color','k'); % Plot nodes
for jj=1:N+1
    distributions={'GLL','lin'};
    dist = distributions{1};
    for j = 1:N+1
        for i = 1:length(z)
            pi(i,j) = Lagrange(j, xi', z(i));
            D(i,j) = derivativeLagrange(j, xi', z(i));
        end
    end
    
    % Plot Lagrange
    plot(z, pi(:,jj)); % Plot cardinal polynomial
    xlim([-1.1, 1.1]); % Set the x-axis limits
    ylim([-0.25 1.1])
    xticks(xi); % Set the x-ticks
    grid on; % Enable the grid
    % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
    xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
    ylabel(['$\ell_i^{(',num2str(N),')}$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
    % axis off
    enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
    xticklabels({'\xi_1','\xi_2','\xi_3','\xi_4','\xi_5','\xi_6','\xi_7','\xi_8' })
    % line([xi(j) xi(j)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line

end
line([[xi(1) xi(end)]],[0 0], 'Color','k', 'LineStyle', '-','LineWidth',1, 'HandleVisibility', 'off')
% legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','$\ell_5$','$\ell_6$','$\ell_7$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
saveas(gcf,['Figures\lagrange_poly_N',num2str(N)],'epsc')
end
set(groot, 'defaultAxesXColor', 'remove', 'defaultAxesYColor', 'remove', 'defaultAxesZColor', 'remove');
    legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','$\ell_5$','$\ell_6$','$\ell_7$','$\ell_8$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','eastoutside','numColumns',1)
fig.Position=[1000 526.6000 341.6000 245.6000]
    saveas(gcf,['Figures\lagrange_poly_legend'],'epsc')



fig=figure;
fig.Position=[746.6000 456.2000 560.0000 267.2000]
hold on
plot_func = 3;
N=3
xi = lglnodes(N);
plot([xi], [zeros(N+1, 1)], 'o', 'Color', defaultColors(1,:)); % Plot nodes
plot([xi], [ones(N+1, 1)], 'o', 'Color', defaultColors(1,:)); % Plot nodes
for jj=1:N+1
    distributions={'GLL','lin'};
    dist = distributions{1};
    for j = 1:N+1
        for i = 1:length(z)
            pi(i,j) = Lagrange(j, xi', z(i));
            D(i,j) = derivativeLagrange(j, xi', z(i));
        end
    end
    
    % Plot Lagrange
    plot(z, pi(:,jj)); % Plot cardinal polynomial
    xlim([-1.1, 1.1]); % Set the x-axis limits
    xticks(xi); % Set the x-ticks
    grid on; % Enable the grid
    % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
    xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
    ylabel(['$\ell_i^{(',num2str(N),')}$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
    % legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
    % axis off
    enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
    xticklabels({'\xi_1','\xi_2','\xi_3','\xi_4' })
    % line([xi(j) xi(j)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line

end
line([[xi(1) xi(end)]],[0 0], 'Color', defaultColors(1,:), 'LineStyle', '-','LineWidth',1)
% legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
saveas(gcf,['Figures\lagrange_poly_N3_appendix'],'epsc')

fig=figure;
fig.Position=[746.6000 456.2000 560.0000 267.2000]
hold on
plot_func = 3;
N=3
xi = lglnodes(N);
plot([xi], [zeros(N+1, 1)], 'o', 'Color', defaultColors(1,:)); % Plot nodes
plot([xi], [zeros(N+1, 1)], 'o', 'Color', defaultColors(1,:)); % Plot nodes
for jj=1:N+1
    distributions={'GLL','lin'};
    dist = distributions{1};
    for j = 1:N+1
        for i = 1:length(z)
            D(i,j) = derivativeLagrange(j, xi', z(i));
        end
    end
    
    % Plot Lagrange
    plot(z, D(:,jj)); % Plot cardinal polynomial
    xlim([-1.1, 1.1]); % Set the x-axis limits
    xticks(xi); % Set the x-ticks
    grid on; % Enable the grid
    % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
    xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
    ylabel(['$\frac{d l_i}{d \xi}$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
    % legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
    % axis off
    enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
    xticklabels({'\xi_1','\xi_2','\xi_3','\xi_4' })
    % line([xi(j) xi(j)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line

end
line([[xi(1) xi(end)]],[0 0], 'Color', defaultColors(1,:), 'LineStyle', '-','LineWidth',1)
% legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
saveas(gcf,['Figures\derivative_lagrange_poly_N3'],'epsc')



fig=figure;
fig.Position=[729 540.2000 581.6000 267.2000]
hold on
plot_func = 3;
N=3;
xi = lglnodes(N);
plot([xi], [zeros(N+1, 1)], 'o', 'Color', defaultColors(1,:)); % Plot nodes
plot([xi], [zeros(N+1, 1)], 'o', 'Color', defaultColors(1,:)); % Plot nodes
for jj=1:N+1
    distributions={'GLL','lin'};
    dist = distributions{1};
    for j = 1:N+1
        for i = 1:length(z)
            DD(i,j) = doublederivativeLagrange(j, xi', z(i));
        end
    end
    
    % Plot Lagrange
    plot(z, DD(:,jj)); % Plot cardinal polynomial
    xlim([-1.1, 1.1]); % Set the x-axis limits
    xticks(xi); % Set the x-ticks
    grid on; % Enable the grid
    % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
    xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
    ylabel(['$\frac{d^2 l_i}{d \xi^2}$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
    % legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
    % axis off
    enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
    xticklabels({'\xi_1','\xi_2','\xi_3','\xi_4' })
    % line([xi(j) xi(j)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line

end
line([[xi(1) xi(end)]],[0 0], 'Color', defaultColors(1,:), 'LineStyle', '-','LineWidth',1)
% legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
saveas(gcf,['Figures\dderivative_lagrange_poly_N3'],'epsc')
%%


[x,w]=lgwt(N-1,-1,1)

fig=figure;
fig.Position=[746.6000 456.2000 560.0000 316]
hold on
plot_func = 3;
N=3;
xi = lglnodes(N);
plot([x], [ones(N-1,1).*2], 'o', 'Color', defaultColors(1,:)); % Plot nodes
plot([x], [ones(N-1, 1).*2], 'o', 'Color', defaultColors(1,:)); % Plot nodes
for jj=1:N+1
    distributions={'GLL','lin'};
    dist = distributions{1};
    for j = 1:N+1
        for i = 1:length(z)
            pi(i,j) = Lagrange(j, xi', z(i));
            D(i,j) = derivativeLagrange(j, xi', z(i));
        end
    end
    
    % Plot Lagrange
    plot(z, pi(:,jj)); % Plot cardinal polynomial
    xlim([-1.1, 1.1]); % Set the x-axis limits
    ylim([-0.25 1.1])
    xticks(x); % Set the x-ticks
    grid on; % Enable the grid
    % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
    ylabel(['$\ell_i^{(3)}$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
    legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
    % axis off
    enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
    xticklabels({'\zeta_1','\zeta_2','\zeta_3','\zeta_4' })
    % line([xi(j) xi(j)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line

end
line([[xi(1) xi(end)]],[0 0], 'Color', defaultColors(1,:), 'LineStyle', '-','LineWidth',1)
legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
saveas(gcf,['Figures\lagrange_poly_N3__zeta'],'epsc')

fig=figure;
fig.Position=[746.6000 456.2000 560.0000 267.2000]
hold on
plot_func = 3;
N=3
xi = lglnodes(N);
plot([x], [ones(N-1, 2)].*20, 'o', 'Color', defaultColors(1,:)); % Plot nodes
% plot([x], [ones(N-1, 2)].*20, 'o', 'Color', defaultColors(1,:)); % Plot nodes
for jj=1:N+1
    distributions={'GLL','lin'};
    dist = distributions{1};
    for j = 1:N+1
        for i = 1:length(z)
            D(i,j) = derivativeLagrange(j, xi', z(i));
        end
    end
    
    % Plot Lagrange
    plot(z, D(:,jj)); % Plot cardinal polynomial
    xlim([-1.1, 1.1]); % Set the x-axis limits
    xticks(x); % Set the x-ticks
    grid on; % Enable the grid
    % ylim([-1.1, 1.1]); % Set the y-axis limits if needed
    xlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', fntsize); % Set the x-label with LaTeX formatting
    ylabel(['$\frac{d l_i}{d \zeta}$'], 'Interpreter', 'latex', 'FontSize', fntsize); % Set the y-label with LaTeX formatting
    % legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
    % axis off
    enhance_plot(0, fntsize, 0, 0, 0); % Enhance plot appearance
    xticklabels({'\zeta_1','\zeta_2','\zeta_3','\zeta_4' })
    % line([xi(j) xi(j)], [0 1],'Color',defaultColors(1,:),'LineWidth',0.7) % Draws a vertical line

end
ylim([-5 5])
line([[xi(1) xi(end)]],[0 0], 'Color', defaultColors(1,:), 'LineStyle', '-','LineWidth',1)
% legend('','','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','', 'Interpreter', 'latex', 'FontSize', fntsize,'location','Northoutside','numColumns',4)
saveas(gcf,['Figures\derivative_lagrange_poly_N3_zeta'],'epsc')




hej=1

function ddl = doublederivativeLagrange(j, xi, x)
% Calculates the second derivative of the j-th Lagrange interpolation polynomial at a point x
%
% Inputs:
%   j - Index of the polynomial
%   xi - Vector of x coordinates
%   x - The point at which to evaluate the derivative
%
% Output:
%   ddl - The value of the second derivative at point x

k = length(xi); % Number of points

% Initialize the second derivative sum
ddl = 0;
for i = 1:k
    if i ~= j
        sum_m = 0;
        term1 = 1 / (xi(j) - xi(i));
        for m = 1:k
            if m ~= j && m ~= i
            term = 1 / (xi(j) - xi(m));
                for n = 1:k
                    if n ~= j && n ~= i && n ~= m
                        term = term * (x - xi(n)) / (xi(j) - xi(n));
                    end
                end
                sum_m = sum_m + term; 
            end
        end
        ddl = ddl + term1 * sum_m;
    end
end
end

function dl = derivativeLagrange(j, xi, x)
% Calculates the j'th Lagrange interpolation polynomial's derivative at a point x
%
% Inputs:
%   j - Derivative of j'th polynomium
%   xi - vector of xi coordinates
%   x - the point at which to evaluate the derivative
%
% Output:
%   dl - the value of the derivative at point z

dl = 0;
k = length(xi);
for i = 1:k
    if i ~= j
        term = 1 / (xi(j) - xi(i));
        for m = 1:k
            if m ~= j && m ~= i
                term = term * (x - xi(m)) / (xi(j) - xi(m));
            end
        end
        dl = dl + term;
    end
end
end

function l = Lagrange(j, xi, x)
% Computes the j-th Lagrange basis polynomial at a point x.
%
% Inputs:
%   j - j'th Lagrange polynomium
%   xi - vector of xi coordinates
%   x - the point at which to evaluate the basis polynomial
%
% Output:
%   l - the value of the i-th basis polynomial at point z

k = length(xi);
l = 1;
for m = 1:k
    if m ~= j
        l = l * (x - xi(m)) / (xi(j) - xi(m));
    end
end
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


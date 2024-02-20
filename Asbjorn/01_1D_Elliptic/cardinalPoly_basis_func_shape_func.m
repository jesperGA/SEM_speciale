% Clean up the workspace and figures to start fresh
clear; close all; clc;

% Add necessary paths for project dependencies
addpath('FEM', 'MESH', 'PLOT', 'TOPOPT', 'Saved_values');

% Initialize variables
N = 4; % Degree of polynomial
[xi, w] = lglnodes(N); % Legendre-Gauss-Lobatto nodes and weights

% Plot Cardinal Polynomials
figure;
for e = 1:N+1
    plot(xi, zeros(N+1, 1), 'o'); % Plot nodes
    hold on;
    xi_lin = linspace(-1, 1, 100); % Fine grid for plotting
    plot(xi_lin, CardinalPolynomial(xi, e, xi_lin), 'Color', 'red'); % Plot cardinal polynomial
    enhance_plot(0,0,0,0,0); % Enhance plot aesthetics
    hold off;
end

% % Plot Legendre Polynomials
% figure;
% plotLegendrePolynomials(N);

% for i = 1 : N+1
%     for j = 1 : N+1
%         d=i==j; % Kronecker delta
%         B(i,j) = B(i,j) + w(j) * J * d;
%         for k = 1 : N+1
%             dli_xi_k = DerivativeLagrangePoly(N, xi, k, i);
%             dlj_xi_k = DerivativeLagrangePoly(N, xi, k, j);
%             A(i,j) = A(i,j) + exp(x(1)+(xi(k)+1)*J) * w(k) * dlj_xi_k * dli_xi_k * 1 / J;
%         end
%     end
% end
% for k = 1:5
%     dli_xi_k = DerivativeLagrangePoly(N, xi, k, 5)
% end
dli_xi_k = DerivativeLagrangePoly(N, 0.5, 1, 5)


function d_ij = d_ij_function(i,j,N,xi)
    d_ij = 0;
    if i==0 && j==0
        d_ij = -1/4*N*(N+1);
    elseif i>=0 && i<=N && j>=0 && j<=N && i~=j
        L_N_xi_i = LegendrePoly(N,xi(i+1));
        L_N_xi_j = LegendrePoly(N,xi(j+1));
        d_ij = L_N_xi_i/L_N_xi_j*(1/(xi(i+1)-xi(j+1)));
    elseif i>=1 && i==j && j<=N-1
        d_ij = 0;
    elseif i==N && j==N
        d_ij = 1/4*N*(N+1);
    end
end

function dl = DerivativeLagrangePoly(N, xi, i, k)
    dl=0;
    for j=1:N+1
        d_ij = d_ij_function(i-1,j-1,N,xi);
        lk = CardinalPolynomial(xi,k,xi(j));
        dl = dl + d_ij*lk;
    end
    
end

% Function to plot Legendre Polynomials
function plotLegendrePolynomials(N)
    for NIdx = 1:5
        xi = linspace(-1, 1, 100); % Fine grid for plotting
        [nodes, ~] = lglnodes(NIdx); % Legendre-Gauss-Lobatto nodes and weights
        plot(nodes, zeros(NIdx+1, 1), 'o'); % Plot nodes
        hold on;
        xticks(nodes); % Set x-ticks at LGL nodes
        grid on; % Enable grid
        L_N = arrayfun(@(xiVal) LegendrePoly(NIdx, xiVal), xi); % Compute Legendre Polynomial values
        plot(xi, L_N); % Plot Legendre Polynomial
        enhance_plot(0,0,0,0,0); % Enhance plot aesthetics
        hold off;
    end
end

% Function to compute Legendre Polynomial of degree N at point xi
function L_N = LegendrePoly(N, xi)
    L = zeros(N+1, 1);
    L(1) = 1;
    L(2) = xi;
    for i = 3:N+1
        n = i - 1;
        L(i) = ((2*n - 1)*xi*L(i - 1) - (n - 1)*L(i - 2)) / n;
    end
    L_N = L(end);
end

% Assume enhance_plot() is defined elsewhere to improve plot aesthetics

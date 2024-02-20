% Clean up the workspace and figures to start fresh
clear; close all; clc;

% Add necessary paths for project dependencies
addpath('FEM', 'MESH', 'PLOT');

% Initialize variables
N = 3; % Degree of polynomial
[xi, w] = lglnodes(N); % Legendre-Gauss-Lobatto nodes and weights
[xi_gl,w_gl]=lgwt(N-2+1,-1,1); 

    % % Compute the corresponding function values
    % y = f(xi);
[xi,w,h]=GetGLL(N+1);
z=linspace(-1,1);
% z=xi
plotto = 2;
for j = 1:N+1
    for i = 1:length(z)
        pi(i,j) = Lagrange(j, xi, z(i));
        D(i,j) = derivativeLagrange(j, xi, z(i));
    end
end

% Plot Lagrange
figure;
for j = 1:plotto
    plot(xi, zeros(N+1, 1), 'o'); % Plot nodes
    hold on;
    plot(z, pi(:,j)); % Plot cardinal polynomial
    enhance_plot(0,0,0,0,0); % Enhance plot aesthetics
    hold off;
end

% Plot Derivative Lagrange
figure;
for j = 1:plotto
    plot(xi, zeros(N+1, 1), 'o'); % Plot nodes
    hold on;
    plot(z, D(:,j)); % Plot cardinal polynomial
    enhance_plot(0,0,0,0,0); % Enhance plot aesthetics
    hold off;
end

z=xi_gl
for j = 1:N+1
    for i = 1:length(z)
        pi_j_zeta_i(i,j) = Lagrange(j, xi, z(i));
        D_j_zeta_i(i,j) = derivativeLagrange(j, xi, z(i));
    end
end

figure(1)
hold on
plot(xi_gl,pi_j_zeta_i(:,end),'o')
enhance_plot(0,0,0,0,0)

figure(2)
hold on
plot(xi_gl,D_j_zeta_i(:,end),'o')
enhance_plot(0,0,0,0,0)

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


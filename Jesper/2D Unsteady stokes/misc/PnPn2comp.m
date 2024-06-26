clear
clc
close all

Pn = readmatrix('PnPn_error.csv');
Pn2 = readmatrix('PnPn-2_error.csv');

DC = get(groot, 'DefaultAxesColorOrder');

figure();

% Create plot objects for each line
semilogy(Pn(:,1), Pn(:,2), '-o','LineWidth',2,'MarkerSize',10,'Color',DC(2,:));
hold on
semilogy(Pn(:,1), Pn(:,3), '-o','LineWidth',2,'MarkerSize',10,'Color',DC(4,:));
semilogy(Pn2(:,1), Pn2(:,2), '-o','LineWidth',2,'MarkerSize',10,'Color',DC(1,:));
semilogy(Pn2(:,1), Pn2(:,3), '-o','LineWidth',2,'MarkerSize',10,'Color',DC(3,:));

set(gca, 'FontSize', 16);

% Adding labels with LaTeX interpreter and specified font sizes
xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| Error \|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 18);

% Adding a legend for each plot line using the plot objects
legend({'$u \in \bf{P}_n$', '$p$ - $P_n$','$u$ - $P_{n-2}$','$p$ - $P_{n-2}$'},...
    'Interpreter', 'latex', 'FontSize', 18,'Location','northoutside','NumColumns',2);
grid on
hold off;
function plotNodalSolution(F,X,U)
    % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    figure;
    scatter3(X(:, 2), X(:, 3), U);
    hold on
    scatter3(X(:,2),X(:,3),F(X(:,2),X(:,3)),'*')
    % Labels and title
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 18);
    zlabel('$u$', 'Interpreter', 'latex', 'FontSize', 18);
    title('Nodal solution', 'Interpreter', 'latex', 'FontSize', 18);
    % axis equal;
    grid on;
    % xlim([-0.1, 1.1]);
    % ylim([-0.1, 1.4]);
    enhance_plot(0, 0, 0, 0, 0);

    % Hold off to stop adding to the current plot
    hold off;
    
    legend('SEM','Analytical')

    % saveas(gcf,'NodalSol','epsc')
end
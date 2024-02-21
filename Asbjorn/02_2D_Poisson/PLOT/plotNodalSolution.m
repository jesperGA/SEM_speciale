function plotNodalSolution(F,X,U)
    figure()
    % Plot convergence
    fntsize = 18;
    defaultColors = get(groot, 'DefaultAxesColorOrder');

    scatter3(X(:, 2), X(:, 3), (U),'LineWidth',2,'MarkerEdgeColor',defaultColors(1,:));

    hold on
    scatter3(X(:,2),X(:,3),(F(X(:,2),X(:,3))),'*','LineWidth',1)
    % Labels and title
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fntsize);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fntsize);
    zlabel('$u$', 'Interpreter', 'latex', 'FontSize', fntsize);
    % axis equal;
    grid on;
    % xlim([-0.1, 1.1]);
    % ylim([-0.1, 1.4]);

    legend('SEM','Analytical','Location','northeast')
    
    enhance_plot(0, fntsize, 20, 0, 2);

    saveas(gcf,'NodalSol','epsc')
end
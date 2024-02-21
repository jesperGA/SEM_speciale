function plotNodalSolution(F,X,U)
    % Plot convergence
    fntsize = 18;
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    figure;
    plot(X(:, 2), U,'o');
    hold on
    plot(linspace(X(1,2),X(end,2)),F(linspace(X(1,2),X(end,2))))
    % Labels and title
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 18);

    % axis equal;
    grid on;
    % xlim([-0.1, 1.1]);
    % ylim([-0.1, 1.4]);
    xticks([0 pi])
    xticklabels({'0','\pi'})

    legend('SEM','Analytical -sin($x$)','Location','north', 'Interpreter', 'latex', 'FontSize', 14)
  
    enhance_plot(0, fntsize, 0, 0, 2);
    
    saveas(gcf,'NodalSol','epsc')
end
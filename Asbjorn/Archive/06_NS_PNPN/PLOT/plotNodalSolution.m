function plotNodalSolution(F,X,U,name,study,n)

        clf
        % Plot convergence
        fntsize = 18;
        defaultColors = get(groot, 'DefaultAxesColorOrder');
        if strcmp(name,'$u_1$') || strcmp(name,'$u_2$')
            scatter3(X(:, 2), X(:, 3), (U),'LineWidth',1,'MarkerEdgeColor',defaultColors(1,:));
        elseif strcmp(name,'$p$')
            scatter3(X(:, 2), X(:, 3), (U),'LineWidth',1,'MarkerEdgeColor',defaultColors(5,:));
        end
        hold on
        scatter3(X(:,2),X(:,3),(F(X(:,2),X(:,3),study.t(n))),'*')
        % Labels and title
        xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fntsize);
        ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fntsize);
        zlabel(name, 'Interpreter', 'latex', 'FontSize', fntsize);
        % title(name, 'Interpreter', 'latex', 'FontSize', 18);
        % axis equal;
        grid on;
        % xlim([-0.1, 1.1]);
        % ylim([-0.1, 1.4]);
        % zlim([-1 1])
        enhance_plot(0, fntsize, 20, 0, 0);

    % saveas(gcf,['Figures\NodalSol_',name],'epsc')
end
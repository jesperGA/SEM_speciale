function plotNodalSolution(F,X,U,name,time)
        clf
        % Plot convergence
        fntsize = 18;
        sizeValue = 100;
        linewidth = 1;
        defaultColors = get(groot, 'DefaultAxesColorOrder');
        if strcmp(name,'$u_1$') || strcmp(name,'$u_2$')
            scatter3(X(:, 2), X(:, 3), (U),'LineWidth',linewidth,'MarkerEdgeColor',defaultColors(1,:),'SizeData', sizeValue);
        elseif strcmp(name,'$p$')
            scatter3(X(:, 2), X(:, 3), (U),'LineWidth',linewidth+1,'MarkerEdgeColor',defaultColors(5,:),'SizeData', sizeValue);
        end
        hold on
        scatter3(X(:,2),X(:,3),(F(X(:,2),X(:,3),time)),'*','SizeData', sizeValue,'LineWidth',linewidth)
        % Labels and title
        xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fntsize);
        ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fntsize);
        zlabel(name, 'Interpreter', 'latex', 'FontSize', fntsize);
        % title(name, 'Interpreter', 'latex', 'FontSize', 18);
        % axis equal;
        grid on;
        % xlim([-0.1, 1.1]);
        % ylim([-0.1, 1.4]);
        % zlim([min(min(U)) max(max(U))])
        % zlim([-0.5 0.15])
        enhance_plot(0, fntsize, 20, 0, 0);

    % saveas(gcf,['Figures\NodalSol_',name],'epsc')
end
function plotConvergence(error,dofs)
    load roenquist_convergence.csv
    
    % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    fig = figure;
    fig.Position = [585 692.2000 525.6000 385.8000];
    semilogy(dofs,error,'o-');
    hold on;
    semilogy(roenquist_convergence(:,1),roenquist_convergence(:,2),'*-');
    
    % Labels and title
    xlabel('$N_t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Error $\left\|u-u_h\right\|_{L^{\infty}, GL}$', 'Interpreter', 'latex', 'FontSize', 18);

    % axis equal;
    grid on;
    xlim([0, 35]);
    ylim([1e-15, 1e0]);
    
    legend('Our Results','R\o nquist ', ...
        'Interpreter', 'latex', 'FontSize', 18, 'Location','southwest' )
    enhance_plot(0, 22, 0, 0, 2);

    saveas(gcf,'convergence','epsc')
end


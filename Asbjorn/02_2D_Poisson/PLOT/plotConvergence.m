function plotConvergence(error)
    load roenquist_convergence_poisson.csv
    
    % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    figure;
    semilogy(5:2:length(error)*2+2,error(2:end),'o-');
    hold on;
    semilogy(roenquist_convergence_poisson(:,1),roenquist_convergence_poisson(:,2),'*-');
    
    % Labels and title
    xlabel('$N_t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Error, $\| u-u_{analytic} \|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 18);
    title('Convergence');
    % axis equal;
    grid on;
    % xlim([-0.1, 1.1]);
    % ylim([-0.1, 1.4]);
    enhance_plot(0, 0, 0, 0, 0);
    
    % Hold off to stop adding to the current plot
    hold off;
    
    % legend('Nodal values E1','Nodal values E2','Nodal values E3','Nodal values E4')
    legend('Our convergence','Rønquist convergence')

    saveas(gcf,'convergence','epsc')
end


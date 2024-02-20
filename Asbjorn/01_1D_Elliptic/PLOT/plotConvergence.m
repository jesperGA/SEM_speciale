function plotConvergence(error,dofs)
    load roenquist_convergence.csv
    
    % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    figure;
    semilogy(dofs,error,'o-');
    hold on;
    semilogy(roenquist_convergence(:,1),roenquist_convergence(:,2),'*-');
    
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


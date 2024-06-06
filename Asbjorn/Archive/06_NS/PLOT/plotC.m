function plotConvergence(mesh,error,Roenquist)
    
    % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
        hold on;
    plot(round(Roenquist(:,1)),Roenquist(:,2),'*');
    plot([0:45],error,'o');
    
    % Labels and title
    xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\lambda_k \frac{2}{\pi N K}$', 'Interpreter', 'latex', 'FontSize', 18);
    % title('Convergence');
    % axis equal;
    grid on;
    % xlim([-0.1, 1.1]);
    % ylim([-0.1, 1.4]);
    
    % legend('Nodal values E1','Nodal values E2','Nodal values E3','Nodal values E4')

    % saveas(gcf,'Figures\convergence','epsc')
end


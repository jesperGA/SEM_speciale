function plotConvergence(mesh,error,Roenquist)
    
    % Plot convergence
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    semilogy(5:2:length(error)*2+2,error(2:end),'o-','Color',defaultColors(5,:));
    hold on;
    semilogy(Roenquist(:,1),Roenquist(:,2),'*-','Color',defaultColors(4,:));
    
    % Labels and title
    xlabel('$N_t$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Error', 'Interpreter', 'latex', 'FontSize', 18);
    % title('Convergence');
    % axis equal;
    grid on;
    % xlim([-0.1, 1.1]);
    % ylim([-0.1, 1.4]);
    
    % legend('Nodal values E1','Nodal values E2','Nodal values E3','Nodal values E4')

    saveas(gcf,'Figures\convergence','epsc')
end


function plotSol2D(mesh,u1,u2,study,t)

    % Plot mesh
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    % figure;
    % plot(  mesh.X(:, 2),   mesh.X(:, 3), '.');
    hold on;
    
    % % Display node numbers
    % for i = 1:length( X)
    %     text(  mesh.X(i, 2),   mesh.X(i, 3), num2str(  mesh.X(i, 1)), ...
    %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    % end
    
    % Labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    % title('Mesh');
    axis equal;
    grid on;
    xticks(unique(sort(  mesh.X(:, 2))));
    yticks(unique(sort(  mesh.X(:, 2))));
    xtickformat('%.2f');
    ytickformat('%.2f');
    constant = 0.5;
    xlim([min(  mesh.X(:,2))-constant, max(  mesh.X(:,2))+constant]);
    ylim([min(  mesh.X(:,3))-constant, max(  mesh.X(:,3))+constant]);
    
    % Plot elements
    for e = 1:size( mesh.IX, 3)
        element_nodes =   mesh.IX(:, :, e);
        x_lower = min(  mesh.X(element_nodes, 2));
        x_upper = max(  mesh.X(element_nodes, 2));
        y_lower = min(  mesh.X(element_nodes, 3));
        y_upper = max(  mesh.X(element_nodes, 3));

        % % Plot box
        % text((x_lower + x_upper) / 2, (y_lower + y_upper) / 2, num2str(e), ...
        %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        %     'Color', 'k', 'FontSize', 12);

        outer_nodes = [element_nodes(1, 1:end), element_nodes(2:end, end)', ...
                       element_nodes(end, end-1:-1:1), element_nodes(end-1:-1:1, 1)'];
        plot(  mesh.X(outer_nodes, 2),   mesh.X(outer_nodes, 3),'--','Color', 'k');
    end

    axis off
    hold on
    % scale_factor = 2; % Adjust the scale factor as needed
    scale_factor = 'off';
    nodes=mesh.bound_u1(:,1);
    nodes=mesh.X(:,1);
    scale = 0.5;
    scale = 1;

    quiver(  mesh.X(nodes,2),   mesh.X(nodes,3), u1(nodes).*scale, u2(nodes).*scale, scale_factor, 'Color', defaultColors(1,:),'LineWidth',1.5)
    quiver(  mesh.X(nodes,2),   mesh.X(nodes,3), study.U1(mesh.X(:,2),mesh.X(:,3),t).*scale, study.U2(mesh.X(:,2),mesh.X(:,3),t).*scale, scale_factor, 'Color', defaultColors(2,:))
    enhance_plot(0, 0, 0.5, 6, 0);
    % saveas(gcf,'Figures\quiver','epsc')
end
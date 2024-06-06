function plotMesh2D_numbering(mesh)

    % Plot mesh
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    fig=figure;
    fig.Position = [651.4000 243.4000 640 606.4000]
    plot(linspace(0,1), ones(1,100),'color', [0.5, 0.5, 0.5]);
    hold on;
    plot([0 0 0.5 0.5 0.5 0.5 1 1 1 0 0], [1 0.5 0.5 1 0 0.5 0.5 1 0 0 0.5], 'color', [0.5, 0.5, 0.5]);
    scatter(mesh.X(:,2),mesh.X(:,3),'.','SizeData', 500, 'MarkerEdgeColor', defaultColors(6,:),'HandleVisibility', 'off')

    fntsize = 22;
    % 
    % % Display node numbers
    % for i = 1:length(mesh.X)
    %     text(mesh.X(i, 2), mesh.X(i, 3), num2str(mesh.X(i, 1)), ...
    %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    % end
    
    % Labels and title
    xlabel('\it x_1');
    ylabel('\it x_2');
    % title('Mesh');
    axis equal;
    grid on;
    xticks(unique(sort(mesh.X(:, 2))));
    yticks(unique(sort(mesh.X(:, 2))));
    xtickformat('%.2f');
    ytickformat('%.2f');
    constant = 0.1;
    xlim([min(mesh.X(:,2))-constant, max(mesh.X(:,2))+constant]);
    ylim([min(mesh.X(:,3))-constant, max(mesh.X(:,3))+constant]);
    
    % Plot elements
    for e = 1:size(mesh.IX, 3)
        element_nodes = mesh.IX(:, :, e);
        x_lower = min(mesh.X(element_nodes, 2));
        x_upper = max(mesh.X(element_nodes, 2));
        y_lower = min(mesh.X(element_nodes, 3));
        y_upper = max(mesh.X(element_nodes, 3));
    
        % Plot box
        text((x_lower + x_upper) / 2 + 0.025, (y_lower + y_upper) / 2+ 0.025, ['\Omega_', num2str(e)], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'k', 'FontSize', 12);
    
        outer_nodes = [element_nodes(1, 1:end), element_nodes(2:end, end)', ...
                       element_nodes(end, end-1:-1:1), element_nodes(end-1:-1:1, 1)'];
        % plot(mesh.X(outer_nodes, 2), mesh.X(outer_nodes, 3), 'k');
        
    end

    % text(1,1.3, '{\it x_2} = 1 + 1/4 sin{\it \pi x_1}', ...
    %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    %         'Color', 'k');
    % text(0, 1, '(0, 1) ', ...
    %         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    %         'Color', 'k');
    % text(1, 0, '(1, 0)', ...
    %         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    %         'Color', 'k');
    text(0, -0.1, '\it{x_1}', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', 'k');
    text(-0.13, 0.05, '\it{x_2}  ', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Color', 'k');
    % Adding an arrow for the x-axis
    % annotation('arrow', [start point x, end point x], [start point y, end point y])
    zero=[0.1 0.1]
    annotation('arrow', [zero(1), zero(1)+0.1], [zero(2), zero(2)], 'Color', 'k');
    annotation('arrow', [zero(1), zero(1)], [zero(2), zero(2)+0.1], 'Color', 'k');
    % zero=[0.67 0.885]
    % annotation('arrow', [zero(1), zero(1)-0.025], [zero(2), zero(2)-0.05], 'Color', 'k');
    box off 
    axis off
    
    
    % Display node numbers
    for i = 1:length(mesh.X)
        text(mesh.X(i, 2), mesh.X(i, 3), num2str(mesh.X(i, 1)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    enhance_plot(0, fntsize, 1, 100, 0);

    saveas(gcf,'2Dmesh_numbering','epsc')

end
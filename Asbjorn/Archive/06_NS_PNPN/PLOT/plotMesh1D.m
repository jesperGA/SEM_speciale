function plotMesh1D(mesh)

    fntsize = 20;
    
    defaultColors = get(groot, 'DefaultAxesColorOrder');

    % Define the nodes (x-coordinates)
    nodes = mesh.X(:,2); % Example nodes
    
    % Define the elements by their connectivity
    elements = mesh.IX(:,2:end-1); % Example elements
    
    % Plotting the mesh
    fig = figure; hold on;
    fig.Position = [744 988.2000 560 120];

    grid on;
    axis off
    
    zero=[0.13 0.58];
    annotation('arrow', [zero(1), zero(1)+0.85], [zero(2), zero(2)], 'Color', 'k');

    % Draw horizontal line for the domain
    xlim([min(nodes), max(nodes)]);
    ylim([-0.1, 0.1]); % Adjust vertical limits for visibility
    line([min(nodes), max(nodes)], [0, 0], 'Color', 'black', 'LineWidth', 1);
    
    % Mark nodes with round markers
    plot(nodes, zeros(size(nodes)), 'o', 'MarkerEdgeColor', defaultColors(1,:), 'MarkerFaceColor', defaultColors(1,:));
    
    for i = 1:size(elements, 1)
        elementNodes = elements(i, :);
        xCoords = nodes(elementNodes);
        
        % Draw thin vertical lines at the start and end of each element
        for j = [1, length(xCoords)]
            line([xCoords(j), xCoords(j)], [-0.05, 0.05], 'Color', 'black', 'LineWidth', 0.5);
        end
    end
    
    height= 0.01;
    % text(pi/4,height, ['\Omega_1'], ...
    % 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    % 'Color', 'k', 'FontSize', fntsize);

    text(-1, -0.05 , '-1 ', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Color', 'k');
    text(1, -0.05, '1', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Color', 'k');
    text(1*1.1, height, '\it{x}', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
            'Color', 'k');
    

    enhance_plot(0,0,0,6,0)
    
    % saveas(gcf,'mesh1D','epsc')
end
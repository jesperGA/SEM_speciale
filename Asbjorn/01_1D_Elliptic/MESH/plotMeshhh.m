function plotMeshhh(mesh)
    % Define the nodes (x-coordinates)
    nodes = mesh.X(:,2); % Example nodes
    
    % Define the elements by their connectivity
    elements = mesh.IX(:,2:end-1); % Example elements
    
    % Plotting the mesh
    fig = figure; hold on;
    fig.Position = [744 988.2000 560 120];

    grid on;
    xticks([0 pi])
    xticklabels({'0','\pi'})
    
    % Remove Y-axis
    set(gca, 'YColor', 'none')
    
    % Draw horizontal line for the domain
    xlim([min(nodes), max(nodes)]);
    ylim([-0.1, 0.1]); % Adjust vertical limits for visibility
    line([min(nodes), max(nodes)], [0, 0], 'Color', 'black', 'LineWidth', 1);
    
    % Mark nodes with round markers
    plot(nodes, zeros(size(nodes)), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');
    
    for i = 1:size(elements, 1)
        elementNodes = elements(i, :);
        xCoords = nodes(elementNodes);
        
        % Draw thin vertical lines at the start and end of each element
        for j = [1, length(xCoords)]
            line([xCoords(j), xCoords(j)], [-0.05, 0.05], 'Color', 'black', 'LineWidth', 0.5);
        end
    end

    hold off;

    saveas(gcf,'mesh1D','epsc')
end

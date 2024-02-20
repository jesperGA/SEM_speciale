function plotMeshh(mesh)
    % Define the nodes (x-coordinates)
    nodes = mesh.X(:,2); % Example nodes
    
    % Define the elements by their connectivity
    elements = mesh.IX(:,2:end-1); % Example elements
    
    % Plotting the mesh
    fig = figure; hold on;
    fig.Position = [744 988.2000 560 120];

    title('Mesh');
    xlabel('x');
    grid on;
    xticks([0 pi])
    xticklabels({'0','\pi'})
    
    % Remove Y-axis
    set(gca, 'YColor', 'none')

    for i = 1:size(elements, 1)
        elementNodes = elements(i, :);
        xCoords = nodes(elementNodes);
        
        % Draw box around element with thick black boundary and no fill
        yOffset = 0.05; % Small vertical offset for visibility
        xBox = [xCoords(1), xCoords(1), xCoords(end), xCoords(end)];
        yBox = [-yOffset, yOffset, yOffset, -yOffset];
        patch(xBox, yBox, 'white', 'EdgeColor', 'black', 'LineWidth', 2);
        
    end

    enhance_plot(0, 0, 0, 0, 0);

    for i = 1:size(elements, 1)
        elementNodes = elements(i, :);
        xCoords = nodes(elementNodes);
        
        % Mark nodes by thin vertical lines
        for j = 1:length(xCoords)
            line([xCoords(j), xCoords(j)], [-yOffset, yOffset], 'Color', 'black', 'LineWidth', 1);
        end
    end

    hold off;

    saveas(gcf,'mesh1D','epsc')
end


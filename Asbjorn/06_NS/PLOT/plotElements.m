function plotMesh2Dgraphic(mesh,ObjectCoords)

    SzData=600;
    lineThk = 3;
    % lineThn = 1.5;



    % Plot mesh
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    fig=figure;
    fig.Position = [97 53.8000 1.0048e+03 778.4000];
    % plot(mesh.X(:, 2), mesh.X(:, 3), '.');
    hold on;
    
    fntsize = 22;
    % Labels and title
    xlabel('\it x_1');
    ylabel('\it x_2');
    axis equal;
    grid on;
    constant = 0.1;
    xlim([min(mesh.X(:,2))-constant, max(mesh.X(:,2))+constant]);
    ylim([min(mesh.X(:,3))-constant, max(mesh.X(:,3))+constant]);
    

    % Initialize arrays for 'p' mesh lines (only interior)
pLinesX = [];
pLinesY = [];

% Initialize arrays to hold line coordinates for efficiency
edgeLinesX = [];
edgeLinesY = [];
interiorLinesX = [];
interiorLinesY = [];

% Process 'v' mesh
for e = 1:size(mesh.IX, 3) % For each element in 'v' mesh

    element_nodes = mesh.IX(:, :, e);
    x_lower = min(mesh.X(element_nodes, 2));
    x_upper = max(mesh.X(element_nodes, 2));
    y_lower = min(mesh.X(element_nodes, 3));
    y_upper = max(mesh.X(element_nodes, 3));

    % Plot box
    % text((x_lower + x_upper) / 2 + 0.025, (y_lower + y_upper) / 2+ 0.025, ['\Omega_{', num2str(e),'}'], ...
    %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    %     'Color', 'k', 'FontSize', 12);

    % Edge lines
    edgeIdx = [1, size(mesh.IX,1); 1, size(mesh.IX,2)]; % Indices for edges (first and last row/column)
    for i = edgeIdx(1,:) % Horizontal edge lines (top and bottom)
        for j = 1:size(mesh.IX,2)-1
            nodeStart = mesh.IX(i, j, e);
            nodeEnd = mesh.IX(i, j+1, e);
            edgeLinesX = [edgeLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
            edgeLinesY = [edgeLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
        end
    end
    for j = edgeIdx(2,:) % Vertical edge lines (left and right)
        for i = 1:size(mesh.IX,1)-1
            nodeStart = mesh.IX(i, j, e);
            nodeEnd = mesh.IX(i+1, j, e);
            edgeLinesX = [edgeLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
            edgeLinesY = [edgeLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
        end
    end

    % % Interior lines
    % for i = 2:size(mesh.IX,1)-1 % Horizontal interior lines
    %     for j = 1:size(mesh.IX,2)-1
    %         nodeStart = mesh.IX(i, j, e);
    %         nodeEnd = mesh.IX(i, j+1, e);
    %         interiorLinesX = [interiorLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
    %         interiorLinesY = [interiorLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
    %     end
    % end
    % for j = 2:size(mesh.IX,2)-1 % Vertical interior lines
    %     for i = 1:size(mesh.IX,1)-1
    %         nodeStart = mesh.IX(i, j, e);
    %         nodeEnd = mesh.IX(i+1, j, e);
    %         interiorLinesX = [interiorLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
    %         interiorLinesY = [interiorLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
    %     end
    % end
end

% Plot lines for 'v' mesh - edges and interior
plot(edgeLinesX, edgeLinesY, 'k-', 'LineWidth', lineThk,'HandleVisibility', 'off'); % Edge lines
% plot(interiorLinesX, interiorLinesY, 'k--', 'LineWidth', lineThn,'HandleVisibility', 'off'); % Interior lines

% % Process 'p' mesh for interior lines
% for e = 1:size(mesh.IXp, 3) % For each element in 'p' mesh
%     for i = 1:size(mesh.IXp,1) % Horizontal lines
%         for j = 1:size(mesh.IXp,2)-1
%             nodeStart = mesh.IXp(i, j, e);
%             nodeEnd = mesh.IXp(i, j+1, e);
%             pLinesX = [pLinesX, [mesh.Xp(nodeStart, 2); mesh.Xp(nodeEnd, 2)]];
%             pLinesY = [pLinesY, [mesh.Xp(nodeStart, 3); mesh.Xp(nodeEnd, 3)]];
%         end
%     end
%     for j = 1:size(mesh.IXp,2) % Vertical lines
%         for i = 1:size(mesh.IXp,1)-1
%             nodeStart = mesh.IXp(i, j, e);
%             nodeEnd = mesh.IXp(i+1, j, e);
%             pLinesX = [pLinesX, [mesh.Xp(nodeStart, 2); mesh.Xp(nodeEnd, 2)]];
%             pLinesY = [pLinesY, [mesh.Xp(nodeStart, 3); mesh.Xp(nodeEnd, 3)]];
%         end
%     end
% end

% % Plot lines for 'p' mesh - interior only, using dotted lines
% plot(pLinesX, pLinesY, 'k:', 'LineWidth',lineThn,'HandleVisibility', 'off');
% 
% scatter(mesh.Xp(:,2),mesh.Xp(:,3),'.', 'SizeData', SzData*1.2, 'MarkerEdgeColor', defaultColors(5,:),'HandleVisibility', 'off')
% scatter(mesh.X(:,2),mesh.X(:,3),'.','SizeData', SzData*1, 'MarkerEdgeColor', defaultColors(1,:),'HandleVisibility', 'off')

% ax = gca; % Get the handle to the current axes
% ax.TickLabelInterpreter = 'latex'; % Set the tick labels to use LaTeX interpreter
% ax.FontSize = 14; % Set the font size for the axes
LX = max(mesh.X(:,2))-min(mesh.X(:,2));
LY = max(mesh.X(:,3))-min(mesh.X(:,3));

if length(ObjectCoords)==4 % square
    x1 = ObjectCoords(1);
    x2 = ObjectCoords(2);
    y1 = ObjectCoords(3);
    y2 = ObjectCoords(4);
    
    % Draw a rectangle on the plot
    rectangle('Position', [x1, y1, x2-x1, y2-y1], 'EdgeColor', defaultColors(2,:), 'LineWidth', 2+lineThk);

elseif length(ObjectCoords)==3 % circle
    % Extract circle parameters
    centerX = ObjectCoords(1);
    centerY = ObjectCoords(2);
    radius = ObjectCoords(3);

    rectangle('Position', [centerX - radius, centerY - radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', defaultColors(2,:), 'LineWidth', lineThk);

    % Calculate the distance from each node to the center of the circle
    distances = sqrt((mesh.X(:,2) - centerX).^2 + (mesh.X(:,3) - centerY).^2);
    
    % Find nodes where the distance is less than or equal to the radius
    Object_nodes = mesh.X(distances <= radius, 1);

    % scatter(mesh.X(Object_nodes,2),mesh.X(Object_nodes,3),'.','SizeData', SzData, 'MarkerEdgeColor', defaultColors(2,:),'HandleVisibility', 'off')

end

% % Add dummy scatter plot for legend
% scatter(NaN, NaN, 'MarkerEdgeColor', defaultColors(1,:), 'LineWidth', 4.5, 'DisplayName', '$u$-grid');
% scatter(NaN, NaN, 'MarkerEdgeColor', defaultColors(5,:), 'LineWidth', 4.5, 'DisplayName', '$p$-grid');
% % Release the hold
% hold off;
% legend('Location','southoutside','Orientation','horizontal','Interpreter','latex',FontSize=fntsize)
hold off;
axis equal;

    % text(1, 1, '(1, 1) ', ...
    %         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    %         'Color', 'k');
    % text(-1, -1, '(-1, -1)', ...
    %         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    %         'Color', 'k');
    % text(1.1, 0, '\it{x_1}', ...
    %         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    %         'Color', 'k');
    % text(0.05, 1.1, '\it{x_2}  ', ...
    %         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    %         'Color', 'k');
    % Adding an arrow for the x-axis
    % annotation('arrow', [start point x, end point x], [start point y, end point y])
    zero=[0.518 0.516];
    arrowlength = 0.39;
    % annotation('arrow', [zero(1), zero(1)+arrowlength], [zero(2), zero(2)], 'Color', 'k');
    % annotation('arrow', [zero(1), zero(1)], [zero(2), zero(2)+arrowlength], 'Color', 'k');
    box off 
    axis off
    
    enhance_plot(0, fntsize, 1, 100, 0);
    % Hold off to stop adding to the current plot
    % hold off;
    K=size(mesh.IX,3);
exportgraphics(gca, ['Figures\meshK',num2str(K),'.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

end
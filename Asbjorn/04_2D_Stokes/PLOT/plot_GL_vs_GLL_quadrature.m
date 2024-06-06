






function plotMesh2Droenquistex(mesh,N)

    % Plot mesh
    defaultColors = get(groot, 'DefaultAxesColorOrder');
    fig=figure;
    fig.Position = [651.4000 419.4000 432 430.4000];
    % plot(mesh.X(:, 2), mesh.X(:, 3), '.');
    hold on;
    
    fntsize = 50;
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
    % 
    % % Plot box
    % text((x_lower + x_upper) / 2 + 0.025, (y_lower + y_upper) / 2+ 0.025, ['\Omega_{', num2str(e),'}'], ...
    %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    %     'Color', 'k', 'FontSize', 12);

    if length(mesh.X)==4
          % Plot box
    text((x_lower + x_upper) / 2 + 0.025, (y_lower + y_upper) / 2+ 0.025, ['\Omega'], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Color', 'k', 'FontSize', 12);
    end

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

    % Interior lines
    for i = 2:size(mesh.IX,1)-1 % Horizontal interior lines
        for j = 1:size(mesh.IX,2)-1
            nodeStart = mesh.IX(i, j, e);
            nodeEnd = mesh.IX(i, j+1, e);
            interiorLinesX = [interiorLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
            interiorLinesY = [interiorLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
        end
    end
    for j = 2:size(mesh.IX,2)-1 % Vertical interior lines
        for i = 1:size(mesh.IX,1)-1
            nodeStart = mesh.IX(i, j, e);
            nodeEnd = mesh.IX(i+1, j, e);
            interiorLinesX = [interiorLinesX, [mesh.X(nodeStart, 2); mesh.X(nodeEnd, 2)]];
            interiorLinesY = [interiorLinesY, [mesh.X(nodeStart, 3); mesh.X(nodeEnd, 3)]];
        end
    end
end

% Plot lines for 'v' mesh - edges and interior
if N==1
    row=1;
elseif N==2
    row=2;
elseif N>2
    row=4;
end
plot(edgeLinesX, edgeLinesY,'Color', defaultColors(row,:),'LineStyle','-' , 'LineWidth',4,'HandleVisibility', 'off'); % Edge lines
plot(interiorLinesX, interiorLinesY,'Color', defaultColors(row,:),'LineStyle','--' , 'LineWidth', 1,'HandleVisibility', 'off'); % Interior lines
% 
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

% Plot lines for 'p' mesh - interior only, using dotted lines
% plot(pLinesX, pLinesY, 'k:', 'LineWidth',1,'HandleVisibility', 'off');

scatter(mesh.X(:,2),mesh.X(:,3),'.','SizeData', 1500, 'MarkerEdgeColor', 'k','HandleVisibility', 'off')
% scatter(mesh.Xp(:,2),mesh.Xp(:,3),'.', 'SizeData', 100, 'MarkerEdgeColor', defaultColors(5,:),'HandleVisibility', 'off')

% ax = gca; % Get the handle to the current axes
% ax.TickLabelInterpreter = 'latex'; % Set the tick labels to use LaTeX interpreter
% ax.FontSize = 14; % Set the font size for the axes

% Add dummy scatter plot for legend
% scatter(NaN, NaN, 'MarkerEdgeColor', defaultColors(1,:), 'LineWidth', 2, 'DisplayName', '$u$-grid');
% scatter(NaN, NaN, 'MarkerEdgeColor', defaultColors(5,:), 'LineWidth', 2, 'DisplayName', '$p$-grid');
% Release the hold
hold off;
% legend('Location','southoutside','Orientation','horizontal','Interpreter','latex',FontSize=fntsize)
hold off;
axis equal;
    % 
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
    % zero=[0.518 0.516];
    % arrowlength = 0.39;
    % annotation('arrow', [zero(1), zero(1)+arrowlength], [zero(2), zero(2)], 'Color', 'k');
    % annotation('arrow', [zero(1), zero(1)], [zero(2), zero(2)+arrowlength], 'Color', 'k');
    box off 
    axis off
    
    enhance_plot(0, fntsize, 2, 100, 0);
    % Hold off to stop adding to the current plot
    % hold off;

end
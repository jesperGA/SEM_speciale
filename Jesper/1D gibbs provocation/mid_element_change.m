clear
close all
clc

GLLS = 4:2:14;

mat = [1,1.2,1.3,1.4;
    2,1.2,1.3,1.4];

L = 1;


subfolder = 'pics_mid_element2';
% Check if the subfolder exists, if not, create it
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

for i =1:numel(GLLS)
    n_GLL = GLLS(i);
    [xi,w,~] = lglnodes(n_GLL-1);

    study.xi = xi;study.w = w;study.n_GLL = n_GLL;
    mesh = [];
    mesh = regular_bragg_grating(2*L,3,n_GLL,mat,study);
    mesh.material = mat(:,1);


    u = zeros(size(mesh.X,1),1);
    u(mesh.X(:,2)>=3*L/4)= 2;
    u(mesh.X(:,2)<3*L/4)= 1;

    [data] = oneD_element_interpolator(mesh,u);
    
    xx = linspace(0,mesh.X(end,2));
    uu = ones(numel(xx),1);
    uu(xx>=3*L/4) = 2;
    
    % xv = mesh.X(:,2);
    
    % plot(xv,u,'ok','LineWidth',2);
    
    figure();
    colors = get(gca,'ColorOrder');
    plot(data(:,1),data(:,2),'LineWidth',2)
    hold on
    s = scatter(mesh.X(:,2),u,100,colors(1,:),'filled');
    plot(xx,uu,'Color',colors(2,:),'LineWidth',2);
    grid on

    % % Define the y-values for the horizontal lines
    y1 = u(1);
    y2 = u(end);

    % Define the x-range for the lines
    x = [0,2*L];

    % Plot the horizontal lines
    % figure();
    % Keep the plot window open to add the second line and texts
    plot(x, [y1 y1], '--k'); % First horizontal line in red
    plot(x, [y2 y2], '--k'); % Second horizontal line in blue
    % hold off;

    fs = 20;

    % Add text annotations to the lines
    % 'HorizontalAlignment' controls the alignment of the text
    text(x(2), y1, ['$\kappa = ' ,num2str(y1),'$'], 'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','FontSize',fs,'Interpreter','latex'); % Text on the right
    text(x(1), y2, ['$\kappa = ' ,num2str(y2),'$'], 'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','FontSize',fs,'Interpreter','latex');  % Text on the left

    title(sprintf('N = %d',n_GLL-1),"FontSize",fs)
    % Set the x-axis limits explicitly
    xlim(x);
    ylim([0,2.5])

    % Customize axis properties
    ax = gca; % Get current axis
    ax.FontSize = fs*0.8; % Set the font size for axis labels

    % Set the x-axis and y-axis ticks
    % x = xlim; % Get the current x-axis limits
    ax.XTick = [x(1), round(mean(x),1), x(2)]; % Set x-axis ticks to start, middle, and end
    % y = ylim; % Get the current y-axis limits
    ax.YTick = [y1, mean([y1,y2]), y2]; % Set y-axis ticks to start, middle, and end
    % Define the subfolder and filename
    plotMeshhh(mesh)

    filename = sprintf('ex_%d.pdf', i);  % Using PNG format, you can change the format as needed
    % Full file path
    fullFilePath = fullfile(subfolder, filename);
    % Save the current figure to the specified path
    % print(gcf, fullFilePath, '-depsc');
    exportgraphics(gca, fullFilePath, 'ContentType', 'vector', 'BackgroundColor', 'none');

    diff2 = max(data(:,2))-y2;
    diff1 = y1-min(data(:,2));

    diffs(:,i) = [diff1;diff2];




end

% figure();
% plot(GLLS,diffs(1,:),'-ok','LineWidth',2)
% xlabel('Number of GLL points in the element: $n_{gll}$','FontSize',fs,'Interpreter','latex')
% ylabel('Overshoot: $\max(p(x))-p_1$','FontSize',fs,'Interpreter','latex')
% ax = gca; % Get current axis
% ax.FontSize = fs*0.8; % Set the font size for axis labels
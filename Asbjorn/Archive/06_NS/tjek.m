dt = 1e-4;
study.t = [0+dt:dt:0.5];
% Example data vectors -- replace these with your actual data
dataSEM = sin(study.t); % Placeholder SEM solution
dataAnalytical = cos(study.t); % Placeholder Analytical solution

figure;
for k = 1:length(study.t)
    clf; % Clear the current figure to redraw
    plot(study.t(1:k), dataSEM(1:k), 'b-', 'LineWidth', 2); % Plot SEM solution up to the k-th point
    hold on; % Hold the plot for additional plotting
    plot(study.t(1:k), dataAnalytical(1:k), 'r--', 'LineWidth', 2); % Plot Analytical solution up to the k-th point
    hold off;
    xlim([0, max(study.t)]); % Keep the x-axis limits static
    % Adjust ylim as necessary
    ylim([min([dataSEM, dataAnalytical]), max([dataSEM, dataAnalytical])]); % Dynamic y-axis limits based on data
    
    % Legend and title to include the current time
    legend('SEM', 'Analytical', 'Location', 'southoutside', 'NumColumns', 2);
    title(sprintf('Time: %0.2f seconds', study.t(k)));
    
    drawnow; % Update the plot
    
    % Capture the plot as an image
    frame = getframe(gcf);
    img = frame2im(frame);
    [imgInd, cmap] = rgb2ind(img, 256); % Convert the image to an indexed image
    
    % Write to the GIF File
    if k == 1
        imwrite(imgInd, cmap, 'timeseries.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imgInd, cmap, 'timeseries.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
function plotStiffnessMatrix(opt)
    % Plot 'stiffnes' matrix
    figure()
    imagesc(opt.A); % This function scales the colors
    colorbar
    
    sum(sum(abs(opt.A)));
    sum(sum(abs(opt.B)));
end
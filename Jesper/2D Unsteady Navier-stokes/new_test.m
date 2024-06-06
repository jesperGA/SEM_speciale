study.T = 8;  % total time in seconds
study.fps = 101;  % frames per second
fdt = 1 / study.fps;  % frame time step
dt = 0.001;  % smaller time step for calculations
tolerance = dt;  % tolerance for frame checking

study.t = 0:dt:study.T;


opt.neqnV = 10;  % Just an example size
opt.neqnP = 5;   % Just an example size
frames = ceil(study.fps * study.T)+1;
opt.U = zeros(opt.neqnV*2, frames);
opt.P = zeros(opt.neqnP, frames);

frame_index = 2;
next_frame_time = fdt;
for i = 2:length(study.t)
    % Simulation logic here
    % e.g., opt.U(:, frame_index) = computeSomething(study.t(i));

    % Check if current time is a frame time
    if study.t(i) >= next_frame_time || i==length(study.t)
        fprintf('Saving frame at t = %f\n', study.t(i));
        % Save data logic
        % opt.U(:, frame_index) = current_data;
        % opt.P(:, frame_index) = current_pressure;
        frame_index = frame_index + 1;
        next_frame_time = next_frame_time + fdt;
    end
end

expected_frames = frames;
if frame_index-1 == expected_frames
    disp('Test Passed: Correct number of frames saved.');
else
    disp(['Test Failed: ', num2str(frame_index-1), ' frames saved, expected ', num2str(expected_frames)]);
end
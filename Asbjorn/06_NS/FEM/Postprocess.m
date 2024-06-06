function opt = Postprocess(mesh, study, opt)
    % Postprocess
    % Checks if data should be saved in NEK5000 format.
%%
    % Time vector
    time = study.t(study.t_steps);

    % Evaluate dissipated energy
    opt.Edis = calculateEnergyDissipation(mesh, study, opt, time);
%%
    vel_mag = sqrt(opt.u1.^2 + opt.u2.^2);

    if strcmp(study.example,'Pipe')
        nodeCoords = [0.615, 0.205];
        % Calculate and display Strouhal number
        opt.St = calculateStrouhal(mesh.X, study, vel_mag, time, nodeCoords);
    
        figure;
        % Plotting velocity vs. time
        plotNodalTimeseries(mesh.X, vel_mag, time, nodeCoords);
        % Plotting pressure vs. time
        plotNodalTimeseries(mesh.Xp, opt.p, time, nodeCoords);
    end

    
    % Save data in NEK5000 format if required
    if study.Save2NEK
        fname = 'nek';
        RefineTimes = 5;
        save2nek5000(fname, opt, mesh, study, RefineTimes);
    end
end

function Edis = calculateEnergyDissipation(mesh, study, opt, time)
    Edis = zeros(opt.nel, length(time));
    [~, w, ~] = GetGLL(study.N + 1);
    for i = 1:length(time)
        for e = 1:opt.nel
            [~, edof_u] = getElementData(mesh.X, mesh.IX, e);
            mu = mesh.Material(1); % Dynamic viscosity
            Edis(e,i) = energyDissipated(mu, opt.u1(edof_u,i), opt.u2(edof_u,i), opt.grad1(:,:,e), opt.grad2(:,:,e), opt.Alpha(edof_u,edof_u), w, opt.J(:,:,e));
        end
    end
    Edis = trapz(time, sum(Edis));
end

% function opt = calculateStrouhal(mesh, study, opt, time)
%     D = 0.1; % Diameter of the obstacle
%     U_m = study.Re / 100; % Mean velocity
%     tol = 1e-12; % Tolerance for selecting specific mesh points
% 
%     % Extract the velocity at the specific point considering the mesh's position
%     u = opt.u1(abs(mesh.X(:,2)-0.615) < tol & abs(mesh.X(:,3)-0.205) < tol, :);
% 
%     % Consider only the latter half of the timeseries for FFT analysis
%     n = length(u);
%     half_n = floor(n / 2);
%     u_half = u(half_n:end);
%     time_half = time(half_n:end);
% 
%     % Perform FFT on the selected half of the velocity data
%     Y = fft(u_half);
%     L = length(u_half);  % Length of the signal
%     P2 = abs(Y/L);  % Two-sided spectrum
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);  % Single-sided spectrum
% 
%     % Frequency domain
%     Fs = 1 / mean(diff(time_half));  % Sampling frequency
%     f = Fs*(0:(L/2))/L;  % Frequency vector
% 
%     % Find the index of the peak frequency
%     [~, I] = max(P1(2:end-1));
%     dominant_frequency = f(I);  % Dominant frequency
% 
%     % Calculate the Strouhal number using the dominant frequency
%     opt.St = D * dominant_frequency / U_m;
% 
%     % Display results
%     disp(['The estimated frequency of oscillation is ', num2str(dominant_frequency), ' Hz']);
%     disp(['The estimated Strouhal number is St = ', num2str(opt.St)]);
% end


function St = calculateStrouhal(X, study, U, time, nodeCoords)
    D = 0.1; % Diameter of the object causing the oscillations
    U_m = study.Re / 100; % Mean velocity based on Reynolds number

    % Calculate the distances from a specific node and find the minimum
    distances = sqrt((X(:,2) - nodeCoords(1)).^2 + (X(:,3) - nodeCoords(2)).^2);
    [~, loc] = min(distances);

    % Extract the velocity time series at the closest node
    u = U(loc, :);

    % Find all peaks in the velocity data
    [pks, locs] = findpeaks(u);
    % 2 peaks pr. period
    pks = pks(2:2:end);
    locs=locs(2:2:end);

    % Using the last three periods
    if length(pks)>4
        last_pks_indices = locs(end-3:end);
    else
        last_pks_indices = 1; % St = NaN
    end
    peak_times = time(last_pks_indices);
    period = diff(peak_times);

    % Calculate the frequency as the inverse of the mean of the three
    % periods
    frequency = mean(1 ./ period);
    St = D * frequency / U_m; % Calculate the Strouhal number

    % Display the calculated frequency and Strouhal number
    disp(['The estimated frequency of oscillation is ', num2str(frequency), ' Hz']);
    disp(['The estimated Strouhal number is St = ', num2str(St)]);
end

function plotNodalTimeseries(X, U, time, nodeCoords)
    distances = sqrt((X(:,2) - nodeCoords(1)).^2 + (X(:,3) - nodeCoords(2)).^2);
    [~, loc] = min(distances);

    u = U(loc, :);
    hold on
    plot(time, u);
    title('Velocity vs. Time');
    xlabel('Time (s)');
    ylabel('Velocity or Pressure (m/s) or (Pa)');
    findpeaks(u, time);
end

function [xy, edof] = getElementData(X, IX, e)
    nen = IX(:,:,e);
    xy = X(nen, 2:3);
    edof = reshape(nen, [], 1);
end

function opt = Postprocess(mesh, study, opt)
    % Postprocess: Performs postprocessing steps including energy dissipation,
    % Strouhal number calculation, and saving results.
    % Inputs:
    %   mesh - Struct containing mesh data
    %   study - Struct containing study parameters
    %   opt - Struct containing solution data
    % Outputs:
    %   opt - Updated struct with postprocessing results

    % Time vector
    time = study.t(study.t_steps);

    % Evaluate dissipated energy
    opt.Edis = calculateEnergyDissipation(mesh, study, opt, time);
    fprintf('The dissipated energy is Phi = %.1f\n', opt.Edis);

    % Evaluate Strouhal number for Duct example
    if strcmp(study.example, 'Duct')
        vel_mag = sqrt(opt.u1.^2 + opt.u2.^2);
        nodeCoords = [0.615, 0.205];
        U_mean = 2 / 3 * max(opt.u1(:, 1)); % Mean velocity based on Reynolds number
        opt.St = calculateStrouhal(mesh.X, study, vel_mag, time, nodeCoords, U_mean);
    end

    % Save data in NEK5000 format if required
    if study.Save2NEK
        fname = 'nek';
        RefineTimes = 5;
        save2nek5000(fname, opt, mesh, study, RefineTimes);
    end

    % Save the mesh, study, and opt structures
    save(fullfile(study.folderName, 'meshStudyOpt.mat'), 'mesh', 'study', 'opt');

    % Create a README file for instructions
    createReadMe(study.folderName);

    fprintf('Study, Opt, and Mesh structures are saved.\nNek folder for plotting in Paraview created.\n');
end

function Edis = calculateEnergyDissipation(mesh, study, opt, time)
    % calculateEnergyDissipation: Calculates the dissipated energy over time
    % Inputs:
    %   mesh - Struct containing mesh data
    %   study - Struct containing study parameters
    %   opt - Struct containing solution data
    %   time - Time vector
    % Outputs:
    %   Edis - Total dissipated energy

    Edis = zeros(opt.nel, length(time));
    [~, w, ~] = GetGLL(study.N + 1);
    for i = 1:length(time)
        for e = 1:opt.nel
            [~, edof_u] = getElementData(mesh.X, mesh.IX, e);
            mu = mesh.Material(1); % Dynamic viscosity
            Edis(e, i) = energyDissipated(mu, opt.u1(edof_u, i), opt.u2(edof_u, i), opt.grad1(:, :, e), opt.grad2(:, :, e), opt.Alpha(edof_u, edof_u), w, opt.J(:, :, e));
        end
    end
    Edis = trapz(time, sum(Edis));
end

function St = calculateStrouhal(X, study, U, time, nodeCoords, U_mean)
    % calculateStrouhal: Calculates the Strouhal number based on velocity data
    % Inputs:
    %   X - Node coordinates
    %   study - Struct containing study parameters
    %   U - Velocity magnitudes
    %   time - Time vector
    %   nodeCoords - Coordinates of the point of interest
    %   U_mean - Mean velocity
    % Outputs:
    %   St - Strouhal number

    D = 2 * study.ObjectCoords(3); % Diameter of the object causing the oscillations

    % Find the closest node to the specified coordinates
    distances = sqrt((X(:, 2) - nodeCoords(1)).^2 + (X(:, 3) - nodeCoords(2)).^2);
    [~, loc] = min(distances);

    % Extract the velocity time series at the closest node
    u = U(loc, :);

    % Find all peaks in the velocity data
    [pks, locs] = findpeaks(u);
    pks = pks(2:2:end);
    locs = locs(2:2:end);

    % Using the last three periods
    if length(pks) > 4
        last_pks_indices = locs(end-3:end);
    else
        last_pks_indices = 1; % St = NaN
    end
    peak_times = time(last_pks_indices);
    period = diff(peak_times);

    % Calculate the frequency as the inverse of the mean of the periods
    frequency = mean(1 ./ period);
    St = D * frequency / U_mean; % Calculate the Strouhal number

    % Display the calculated frequency and Strouhal number
    fprintf('The estimated frequency of oscillation is %.2f Hz\n', frequency);
    fprintf('The estimated Strouhal number is St = %.2f\n', St);
end

function [xy, edof] = getElementData(X, IX, e)
    % getElementData: Extracts element data for the given element index
    % Inputs:
    %   X - Node coordinates
    %   IX - Element connectivity matrix
    %   e - Element index
    % Outputs:
    %   xy - Coordinates of the element nodes
    %   edof - Element degrees of freedom

    nen = IX(:, :, e);
    xy = X(nen, 2:3);
    edof = reshape(nen, [], 1);
end

function createReadMe(folderName)
    % createReadMe: Creates a README file with instructions for visualization
    % Inputs:
    %   folderName - Name of the folder where the README file will be created

    fileName = fullfile(folderName, 'README.txt');
    message = 'To get solution into Paraview, simply drag the nek\\nek.nek5000 file into Paraview.';

    % Open the file for writing
    fileID = fopen(fileName, 'w');

    % Write the message to the file
    fprintf(fileID, '%s\n', message);

    % Close the file
    fclose(fileID);
end

function opt = Controller(mesh, study)
    % Controller manages the flow of finite element analysis and topology optimization
    % by coordinating the assembly, solving, and postprocessing steps.
    % Inputs:
    %   mesh - Struct containing the mesh details
    %   study - Struct containing the study parameters
    % Outputs:
    %   opt - Struct containing optimization results and performance metrics

    % Display problem size information
    ndofs = 2 * size(mesh.X, 1) + size(mesh.Xp, 1);
    fprintf('Number of elements: %d\n', size(mesh.IX, 3));
    fprintf('Number of degrees of freedom (DOFs): %d\n', ndofs);

    % Initialize the output structure
    opt = [];

    % Calculate and display Reynolds number for specific examples
    if strcmp(study.example, 'Duct') || strcmp(study.example, 'LidDriven')
        opt.Re = calculateReynolds(mesh, study);
        fprintf('Reynolds number is Re = %d\n\n', round(opt.Re));
    end

    % Assembly phase
    tic;  % Start timing assembly
    disp('Starting Assembly...');
    opt = Assembly(mesh, study, opt);
    opt.assemblyTime = toc;  % Stop timing assembly
    fprintf('Assembly completed in: %.1f seconds\n\n', opt.assemblyTime);

    % Solution phase
    ticSol = tic;  % Start timing solution
    disp('Starting Solver...');
    opt = Solver(mesh, study, opt);
    opt.solutionTime = toc(ticSol);  % Stop timing solution
    fprintf('\nSolution completed in: %.1f seconds\n\n', opt.solutionTime);

    % Post-processing phase
    tic;  % Start timing post-processing
    disp('Starting Postprocess...');
    opt = Postprocess(mesh, study, opt);
    opt.postProcTime = toc;  % Stop timing post-processing
    fprintf('Postprocessing completed in: %.1f seconds\n\n', opt.postProcTime);
end

function Re = calculateReynolds(mesh, study)
    % calculateReynolds Computes the Reynolds number based on the study parameters
    % Inputs:
    %   mesh - Struct containing the mesh details
    %   study - Struct containing the study parameters
    % Outputs:
    %   Re - Calculated Reynolds number

    mu = mesh.Material(1);
    rho = mesh.Material(2);

    if strcmp(study.example, 'Duct')
        L = 2 * study.ObjectCoords(3);
        u_mean = (2/3) * max(study.U1(mesh.X(:, 2), mesh.X(:, 3)));
        Re = (rho * L * u_mean) / mu;
    elseif strcmp(study.example, 'LidDriven')
        L = max(mesh.X(:, 2)) - min(mesh.X(:, 2));
        u = max(study.U1(mesh.X(:, 2), mesh.X(:, 3)));
        Re = (rho * L * u) / mu;
    end
end

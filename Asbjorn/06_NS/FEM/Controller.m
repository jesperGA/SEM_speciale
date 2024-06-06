function opt = Controller(mesh, study)
    % Controls the flow of finite element analysis and topology optimization
    % by coordinating the assembly, solving, and postprocessing steps.

    % Display problem size information
    ndofs = 2 * size(mesh.X, 1) + size(mesh.Xp, 1);
    fprintf('Number of elements: %d\n', size(mesh.IX, 3));
    fprintf('Number of degrees of freedom (DOFs): %d\n', ndofs);

    % Initialize the output structure
    opt = [];

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
    opt.postProcTime = toc;  % Stop timing post-processings
    fprintf('Postprocessing completed in: %.1f seconds\n\n', opt.postProcTime);
end

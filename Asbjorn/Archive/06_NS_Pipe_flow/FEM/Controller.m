function opt = Controller(mesh,study)
% Main driver for fea and topopt
% opt = Controller(mesh,study)
%
% Input containers:
% mesh.X            = nodal coordinates
% mesh.IX           = topology table
% mesh.Material:    = list of materials (e,[thk E nu dens])
% mesh.mu:          = List of elemental design variables
%
% study.maxit       = Max number of optimization iterations
% study.p           = Penalization factor
% study.N           = Element polynomial order
% study.t_int_type  = Type of time integration
% study.n_t         = Number of time steps
% study.t           = List of time steps
% study.pulse       = f(t) - value of Gaussian pulse for each time step
% study.optimize    = 0, 1
%
% Output container
% opt.mu, P, K, M, D, u, ud, udd, and more... check it out !

% Problem size
fprintf('Number of elements: %d\n', size(mesh.IX,3));
fprintf('Number of dofs: %d\n', 2*size(mesh.X,1)+size(mesh.Xp,1));

%% Initialize empty output
opt = [];

tic
% Assembly
gcp;
disp('ASSEMBLY')
opt = Assembly(mesh, study, opt);
assemblyTime = toc;  % End timer
% Display assembly time
fprintf('Time for assembly: %.1f seconds\n', assemblyTime); fprintf('\n');

tic
% Solve the system
opt = Solver(mesh, study, opt);
solutionTime = toc;    
fprintf('\nTime for solution: %.1f seconds\n', solutionTime);fprintf('\n')

tic
% Postproc for Element sol
opt = Postprocess(mesh,study,opt);
postProcTime = toc;    
fprintf('POSTPROCESS \nTime for Postprocess: %.1f seconds\n', postProcTime);fprintf('\n')


end

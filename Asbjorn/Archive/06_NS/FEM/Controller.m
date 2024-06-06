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
fprintf('Number of elements: %i \nNumber of dofs: %i \n', size(mesh.IX,3),size(mesh.X,1))

%% Initialize empty output
opt = [];

% Perform assembly
opt = Assembly(mesh, study, opt);

tic
% Solve the system
opt = Solver(mesh, study, opt);
toc

%% Postproc for energies, stresses, strains
% opt = Postprocess(mesh,study,opt);

end

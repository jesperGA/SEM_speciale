% Define the size of the matrix
n = 10000; % Size of the matrix

% Generate random diagonal matrix
A_diag = diag(rand(n, 1));

% Generate random symmetric matrix
A_sym = rand(n, n);
A_sym = 0.5*(A_sym + A_sym'); % Make it symmetric

% Generate random vector f
f = rand(n, 1);

% Solve the system Au = f for diagonal matrix
tic; % Start timer
u_diag = A_diag \ f;
time_diag = toc; % End timer

% Solve the system Au = f for symmetric matrix
tic; % Start timer
u_sym = A_sym \ f;
time_sym = toc; % End timer

% Display solution times
fprintf('Time taken to solve for diagonal matrix: %.6f seconds\n', time_diag);
fprintf('Time taken to solve for symmetric matrix: %.6f seconds\n', time_sym);

% Check the correctness of solutions (optional)
fprintf('\nChecking correctness of solutions...\n');
fprintf('Error for diagonal matrix: %.6e\n', norm(A_diag * u_diag - f));
fprintf('Error for symmetric matrix: %.6e\n', norm(A_sym * u_sym - f));
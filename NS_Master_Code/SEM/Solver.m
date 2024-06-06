function opt = Solver(mesh, study, opt)
    % Solver: Computes the solution for a given SEM problem setup.
    % Inputs:
    %   mesh - Struct containing mesh data
    %   study - Struct containing study parameters
    %   opt - Struct for storing assembled matrices and vectors
    % Outputs:
    %   opt.u1, opt.u2, opt.p - Solutions for velocities and pressure

    % Time step
    dt = study.t(2) - study.t(1);
    % Material properties
    mu = mesh.Material(1);
    rho = mesh.Material(2);

    % Define the system matrices based on the time integration method
    if strcmp(study.timeInt, 'FullyExplicit_FD1_AB3')
        H1 = opt.LHS_BC(1:opt.neqn_u, 1:opt.neqn_u);
        H2 = opt.LHS_BC(opt.neqn_u+1:2*opt.neqn_u, opt.neqn_u+1:2*opt.neqn_u);
        D1 = -opt.LHS_BC(2*opt.neqn_u+1:end, 1:opt.neqn_u);
        D2 = -opt.LHS_BC(2*opt.neqn_u+1:end, opt.neqn_u+1:2*opt.neqn_u);
        E = D1 * (H1 \ D1') + D2 * (H2 \ D2');
    end

    % Preallocate arrays for solution
    U_out = zeros(opt.neqn, length(study.t_steps));
    g_u1 = sparse(opt.neqn_u, 1);
    g_u2 = sparse(opt.neqn_u, 1);
    g_p = sparse(opt.neqn_p, 1);

    % Get coefficients for time integration based on specified method
    if strcmp(study.timeInt, 'BDF1_AB3') || strcmp(study.timeInt, 'FullyExplicit_FD1_AB3')
        beta = [1, 0, 0];
        alpha = [23/12, -4/3, 5/12];
    elseif strcmp(study.timeInt, 'BDF3_AB3')
        beta = [3, -3/2, 1/3];
        alpha = [3, -3, 1];
    end
    
    % Derivative Lagrange
    [~, ~, D_hat] = GetGLL(study.N + 1);

    % Initialize velocities for u1, u2
    u1 = sparse(study.U1(mesh.X(:, 2), mesh.X(:, 3), 0));
    u2 = sparse(study.U2(mesh.X(:, 2), mesh.X(:, 3), 0));
    u1_old = sparse(study.U1(mesh.X(:, 2), mesh.X(:, 3), 0 - dt));
    u2_old = sparse(study.U2(mesh.X(:, 2), mesh.X(:, 3), 0 - dt));
    u1_oldold = sparse(study.U1(mesh.X(:, 2), mesh.X(:, 3), 0 - 2 * dt));
    u2_oldold = sparse(study.U2(mesh.X(:, 2), mesh.X(:, 3), 0 - 2 * dt));
    U_out(1:opt.neqn_u, 1) = u1;
    U_out(opt.neqn_u+1:2*opt.neqn_u, 1) = u2;

    % Initialize convection matrices
    ldof = (study.N + 1) ^ 2;
    [BD1, BD2] = precomputeGlobalCmatr(opt.nel, opt.b, D_hat, mesh.L1, mesh.L2);
    [Cu1, Cu2] = assembleC(full(u1), full(u2), BD1, BD2, mesh.IX, opt.neqn_u, opt.nel);
    [Cu1old, Cu2old] = assembleC(full(u1_old), full(u2_old), BD1, BD2, mesh.IX, opt.neqn_u, opt.nel);
    [Cu1oldold, Cu2oldold] = assembleC(full(u1_oldold), full(u2_oldold), BD1, BD2, mesh.IX, opt.neqn_u, opt.nel);

    % Precompute system matrix decomposition based on the time integration method
    if strcmp(study.timeInt, 'BDF1_AB3') || strcmp(study.timeInt, 'BDF3_AB3')
        Dsysmat = decomposition(opt.LHS_BC);
    elseif strcmp(study.timeInt, 'FullyExplicit_FD1_AB3')
        if study.N > 6
            Dsysmat = decomposition(-E, 'banded');
        else
            Dsysmat = decomposition(-E);
        end
    end

    % Check if boundary conditions are time-independent
    time_independent = length(strfind(func2str(study.U1), 't')) <= 1 && length(strfind(func2str(study.U2), 't')) <= 1;
    boundaryNodes = mesh.bound_u1(:, 1);
    free_dofs = diag(opt.Null);

    if time_independent
        g_u1(boundaryNodes) = study.U1(mesh.X(boundaryNodes, 2), mesh.X(boundaryNodes, 3), 0);
        g_u2(boundaryNodes) = study.U2(mesh.X(boundaryNodes, 2), mesh.X(boundaryNodes, 3), 0);
        BC = [g_u1; g_u2; g_p];
        BCindex = find(BC);
    end

    prev_mes=0;
    % Solver loop for each time step
    for n = 2:length(study.t)
        % Construct right-hand side (RHS) for momentum equations
        if strcmp(study.timeInt, 'BDF1_AB3') || strcmp(study.timeInt, 'BDF3_AB3')
            g1 = opt.B * (opt.f_u1 + rho / dt * (beta(1) * u1 + beta(2) * u1_old + beta(3) * u1_oldold)) - rho * (alpha(1) * Cu1 + alpha(2) * Cu1old + alpha(3) * Cu1oldold);
            g2 = opt.B * (opt.f_u2 + rho / dt * (beta(1) * u2 + beta(2) * u2_old + beta(3) * u2_oldold)) - rho * (alpha(1) * Cu2 + alpha(2) * Cu2old + alpha(3) * Cu2oldold);

            % Apply boundary conditions and adjust RHS
            if ~time_independent
                g_u1(boundaryNodes) = study.U1(mesh.X(boundaryNodes, 2), mesh.X(boundaryNodes, 3), study.t(n));
                g_u2(boundaryNodes) = study.U2(mesh.X(boundaryNodes, 2), mesh.X(boundaryNodes, 3), study.t(n));
                BC = [g_u1; g_u2; g_p];
                BCindex = find(BC);
            end

            RHS = [g1; g2; sparse(opt.neqn_p, 1)] - opt.LHS * BC;
            RHS(~free_dofs) = 0;
            RHS(BCindex) = BC(BCindex); % Enforce boundary conditions

            % Solve the system
            U = Dsysmat \ RHS;

            % Update for next time step
            u1_oldold = u1_old;
            u1_old = u1;
            u1 = U(1:opt.neqn_u);
            u2_oldold = u2_old;
            u2_old = u2;
            u2 = U(opt.neqn_u+1:2*opt.neqn_u);
            p = U(2*opt.neqn_u+1:end);

        elseif strcmp(study.timeInt, 'FullyExplicit_FD1_AB3')
            g1 = -mu * opt.A * u1 + opt.B * (opt.f_u1 + rho / dt * (beta(1) * u1)) - rho * (alpha(1) * Cu1 + alpha(2) * Cu1old + alpha(3) * Cu1oldold);
            g2 = -mu * opt.A * u2 + opt.B * (opt.f_u2 + rho / dt * (beta(1) * u2)) - rho * (alpha(1) * Cu2 + alpha(2) * Cu2old + alpha(3) * Cu2oldold);

            % Apply boundary conditions and adjust RHS
            if ~time_independent
                g_u1(boundaryNodes) = study.U1(mesh.X(boundaryNodes, 2), mesh.X(boundaryNodes, 3), study.t(n));
                g_u2(boundaryNodes) = study.U2(mesh.X(boundaryNodes, 2), mesh.X(boundaryNodes, 3), study.t(n));
                BC = [g_u1; g_u2; g_p];
                BCindex = find(BC);
            end

            RHS = [g1; g2; zeros(opt.neqn_p, 1)] - opt.LHS * BC;
            RHS(~free_dofs) = 0;
            RHS(BCindex) = BC(BCindex); % Enforce boundary conditions

            RHS1 = sparse(RHS(1:opt.neqn_u));
            RHS2 = sparse(RHS(opt.neqn_u+1:2*opt.neqn_u));
            RHS3 = sparse(RHS(2*opt.neqn_u+1:end));

            p = Dsysmat \ (D1 * (H1 \ RHS1) + D2 * (H2 \ RHS2) + RHS3);
            u1 = H1 \ (RHS1 + D1' * p);
            u2 = H2 \ (RHS2 + D2' * p);
        end

        % Update convection matrices for next time step
        Cu1oldold = Cu1old;
        Cu2oldold = Cu2old;
        Cu1old = Cu1;
        Cu2old = Cu2;
        [Cu1, Cu2] = assembleC(full(u1), full(u2), BD1, BD2, mesh.IX, opt.neqn_u, opt.nel);

        % Store solution
        if any(study.t(n) == study.t(study.t_steps))
            i = find(study.t(n) == study.t(study.t_steps));
            U_out(:, i) = [u1; u2; p];
        end

        if mod(n, 1000) == 0
            [prev_mes] = statusMsg(n, study.t, prev_mes);
        end

        if any(isnan([u1]))
            disp('The vector contains NaN values.');
            [opt.u1, opt.u2, opt.p] = output(U_out, opt.neqn_u, mesh.Xp, study.P, study.t);
            return;
        end
    end

    [opt.u1, opt.u2, opt.p] = output(U_out, opt.neqn_u, mesh.Xp, study.P, study.t(study.t_steps));
end

function [u1, u2, p] = output(U_out, neqn_u, Xp, P, time)
    % Format the output
    u1 = U_out(1:neqn_u, :);
    u2 = U_out(neqn_u+1:2*neqn_u, :);
    p = U_out(2*neqn_u+1:end, :);
    p = p - mean(p - P(Xp(:, 2), Xp(:, 3), time));
end

function [prev_mes] = statusMsg(n, t, prev_mes)
    % Display the progress status message
    status_percent = (n / (length(t) - 1)) * 100;
    mess = sprintf('%3.2f percent', status_percent);
    back_space = repmat('\b', 1, prev_mes);
    fprintf(back_space);
    fprintf(mess);
    prev_mes = length(mess);
end

function [BD1, BD2] = precomputeGlobalCmatr(numElems, Be, D_hat, L1, L2)
    % Precompute indices and coefficients for constructing global convection matrices.
    % Inputs:
    %   numElems - Number of elements
    %   Be - Element matrices
    %   D_hat - Derivative matrix
    %   L1, L2 - Element length scales
    % Outputs:
    %   BD1, BD2 - Precomputed global convection matrices

    Imat = eye(length(D_hat));
    D1 = kron(Imat, D_hat) * 2;
    D2 = kron(D_hat, Imat) * 2;

    BD1 = zeros(size(Be));
    BD2 = zeros(size(Be));

    for e = 1:numElems
        % Compute scales for element-level convection matrix based on constant parts
        BD1(:, :, e) = Be(:, :, e) * (D1 / L1(e));
        BD2(:, :, e) = Be(:, :, e) * (D2 / L2(e));
    end
end

function [Cu1, Cu2] = assembleC(u1, u2, BD1, BD2, IX, neqn_u, nel)
    % Assemble global convection matrices
    % Inputs:
    %   u1, u2 - Velocity components
    %   BD1, BD2 - Precomputed matrices
    %   IX - Element connectivity
    %   neqn_u - Number of equations for velocity
    %   nel - Number of elements
    % Outputs:
    %   Cu1, Cu2 - Global convection matrices

    Cu1 = zeros(neqn_u, 1);
    Cu2 = zeros(neqn_u, 1);

    for e = 1:nel
        nen = IX(:, :, e);
        nen = nen(:);

        % Accumulate the contributions into the global result vectors
        Cu1(nen) = Cu1(nen) + (u1(nen) .* BD1(:, :, e) + u2(nen) .* BD2(:, :, e)) * u1(nen);
        Cu2(nen) = Cu2(nen) + (u1(nen) .* BD1(:, :, e) + u2(nen) .* BD2(:, :, e)) * u2(nen);
    end
end

% function [Cu1, Cu2] = assembleC(u1, u2, BD1, BD2, IX, neqn_u, nel)
% 
%     % PARFOR NOT FASTER (at least for Coarse meshes)
% 
% % Initialize local accumulators for each worker
% Cu1 = zeros(neqn_u, 1);
% Cu2 = zeros(neqn_u, 1);
% 
% parfor e = 1:nel
%     nen = IX(:,:,e);
%     nen = nen(:);
% 
%     % Local accumulators for this iteration
%     localCu1 = zeros(neqn_u, 1);
%     localCu2 = zeros(neqn_u, 1);
% 
%     % Accumulate the contributions into the local result vectors
%     localCu1(nen) = (u1(nen) .* BD1(:,:,e) + u2(nen) .* BD2(:,:,e)) * u1(nen);
%     localCu2(nen) = (u1(nen) .* BD1(:,:,e) + u2(nen) .* BD2(:,:,e)) * u2(nen);
% 
%     % Aggregate the local results to the global accumulators
%     Cu1 = Cu1 + localCu1;
%     Cu2 = Cu2 + localCu2;
% end
% 
% end


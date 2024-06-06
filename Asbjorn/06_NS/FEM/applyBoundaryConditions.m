function [LHS, RHS] = applyBoundaryConditions(LHS, RHS, Null, g, neqn)
    % applyBoundaryConditions: Applies boundary conditions to the system matrices
    % 
    % Input:
    % opt - Struct containing system matrices and boundary conditions
    % 
    % Output:
    % opt - Updated struct with modified system matrices and load vector

    % Identify free degrees of freedom
    free_dofs = diag(Null);

    % Create a diagonal sparse matrix of boundary values
    G = spdiags(g, 0, neqn, neqn);

    % Modify the load vector
    RHS = RHS - sum(LHS * G, 2);       % Adjust for stiffness matrix
    RHS(~free_dofs) = 0;
    RHS(find(g)) = g(find(g));   % Enforce boundary values on the load vector

    % Modify the stiffness matrix for boundary conditions
    LHS = Null' * LHS * Null - (Null - speye(neqn, neqn));

end
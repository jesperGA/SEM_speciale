function [p, convergenceInfo] = pcgUzawa(S, Bf1, Bf2, Fp, B_p, p0, A1, A2, D1, D2, tol, max_iter)
%pcgUzawa Solves for the pressure p using a preconditioned conjugate gradient method.
% 
% Syntax:
% [p, convergenceInfo] = pcgUzawa(S, Bf1, Bf2, Fp, B_p, p0, A1, A2, D1, D2, tol, max_iter)
%
% Inputs:
% S, Bf1, Bf2, Fp, A1, A2, D1, D2 - Matrix operators
% p0 - Initial guess for the pressure
% B_p - Precondition matrix
% tol - Tolerance for the convergence
% max_iter - Maximum number of iterations
%
% Outputs:
% p - Computed pressure vector
% convergenceInfo - Struct containing convergence status and iteration details

% Initialize variables
u = zeros(size(p0));
RHS = -D1 * (A1 \ Bf1) - D2 * (A2 \ Bf2) - Fp;
r_old = -RHS + S * p0;
q_old = B_p \ r_old;
p = q_old;

% Main PCG iteration loop
for m = 1:max_iter
    y1 = D1'*p;
    y2 = D2'*p;

    z1 = A1\y1; 
    z2 = A2\y2;

    Sp = D1*z1 + D2*z2;

    a = -(dot(q_old, r_old) / dot(p,Sp));
    u = u + a * p;
    r = r_old + a * Sp;
    q = B_p \ r;
    b = dot(q, r) / dot(q_old, r_old);
    p = q + b * p;

    ch = norm(r);

    r_old = r; 
    q_old = q;

    if m>=max_iter || ch < tol
        p = u;
        if m>=max_iter
            warning('pcg not converged succesfully. Max iter reached')
        end
        convergenceInfo = itermsg('pcg', tol, max_iter, m, 0, m, ch);
        break
    end
end

end

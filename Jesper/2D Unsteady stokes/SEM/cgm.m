function [x,stat]=cgm(A,b,x0)
% The conjugate gradient method to solve Ax=b or min 1/2x^TAx-b^Tx
%
% Inputs:
% A: the system matrix of the linear system, or the Hessian of the
% minimization problem.
% b: the right hand side of the linear system, or the vector in
% linear term for the minimization problem.
% x0: the starting point.
%
% Outputs:
% x: the solution, if the algorithm converges.
% stat: a structure array.
% stat.converged shows if the algorithm converges. 1: converged; 0: not.
% stat.iter gives the number of iterations.
% stat.X saves all iterate including the starting point.
% stat.resd saves the 2-norm of all residules, i.e. \|b-Ax\|_2.

% Solver settings and info
maxit = 100;
tol = 5.0e-6;

stat.converged = false; % converged
stat.iter = 0; % number of iterations

% Initial iteration
x = x0;
it = 0;
r = b-A*x;
p = r;
norm_r = norm(r);
converged = false;

% Store data for plotting
stat.X = x;
stat.resd = norm_r; % norm of residuals

% Main loop of conjugate gradient
while ~converged && (it < maxit)
    it = it+1;
    
    Ap = A*p;
    alpha = norm_r^2/(p'*Ap);
    x = x+alpha*p;
    r = r-alpha*Ap;
    norm_r_New = norm(r);
    beta = (norm_r_New/norm_r)^2;
    p = r+beta*p;
    
    norm_r = norm_r_New;    
    converged = (norm_r <= tol);
    
    % Store data for plotting
    stat.X  = [stat.X  x];
    stat.resd  = [stat.resd norm_r];
end

% Prepare return data
if ~converged
    x = [];
end
stat.converged = converged;
stat.iter = it;
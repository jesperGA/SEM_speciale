function opt = Solver(mesh,study,opt)
% function opt = Solver(mesh,study,p)
%
% Input
% study.analysis =  'static'
%                   ! More will come from next week !
%
% study.X with X = neig (int), omega (float(:))
%
% mesh: mesh data (X,IX,bound,....)
% opt: Matrices (K,C,M,P,Null)
%
% Output
% opt.U , opt.w , (opt.pdof)

A1 = opt.LHS(1:opt.neqn_u,1:opt.neqn_u);
A2 = opt.LHS(opt.neqn_u+1:2*opt.neqn_u,opt.neqn_u+1:2*opt.neqn_u);
D1 = -opt.LHS(2*opt.neqn_u+1:end,1:opt.neqn_u);
D2 = -opt.LHS(2*opt.neqn_u+1:end,opt.neqn_u+1:2*opt.neqn_u);

RHS1 = opt.RHS(1:opt.neqn_u);
RHS2 = opt.RHS(opt.neqn_u+1:2*opt.neqn_u);
RHS3 = opt.RHS(2*opt.neqn_u+1:end);

S = D1 * (A1 \ D1') + D2 * (A2 \ D2');
RHS = - D1 * (A1 \ RHS1) - D2 * (A2 \ RHS2) - RHS3;

S = S/2 + S'/2;

p0 = zeros(opt.neqn_p,1);

if strcmp(study.solver,'Direct')==1

    tic;
    U = opt.LHS \ (opt.RHS);
    time = toc;
    fprintf('Solution found in %f s\n',time);

    opt.u1 = U(1:opt.neqn_u);
    opt.u2 = U(opt.neqn_u+1:2*opt.neqn_u);
    opt.p = U(2*opt.neqn_u+1:end);

elseif strcmp(study.solver,'pcg')==1
    
    tol = 1e-10;
    opt.p = pcg(S, RHS, tol, [], opt.M, [] , []);
    opt.u1 = A1 \ (D1' * opt.p + RHS1);
    opt.u2 = A2 \ (D2' * opt.p + RHS2);

elseif strcmp(study.solver,'Uzawa')==1
    
    [opt.p] = pcgUzawa(S, RHS, opt.M, p0, A1, A2, D1, D2);
    opt.u1 = A1 \ (D1' * opt.p + RHS1);
    opt.u2 = A2 \ (D2' * opt.p + RHS2);

elseif strcmp(study.solver,'UzawamodJGA')==1

    [opt.p] = pcg_mod(S, RHS, opt.M, p0, A1, A2, D1, D2);
    opt.u1 = A1 \ (D1' * opt.p + RHS1);
    opt.u2 = A2 \ (D2' * opt.p + RHS2);

end

sum(abs(A1*opt.u1-D1'*opt.p-RHS1))
sum(abs(A2*opt.u2-D2'*opt.p-RHS2))
sum(abs(-D1*opt.u1-D2*opt.u2-RHS3))

end
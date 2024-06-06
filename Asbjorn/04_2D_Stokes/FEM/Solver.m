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

Bf1 = opt.RHS(1:opt.neqn_u);
Bf2 = opt.RHS(opt.neqn_u+1:2*opt.neqn_u);
Fp = opt.RHS(2*opt.neqn_u+1:end);

S = D1 * (A1 \ D1') + D2 * (A2 \ D2');
RHS = - D1 * (A1 \ Bf1) - D2 * (A2 \ Bf2) - Fp;

S = S/2 + S'/2;

p0 = zeros(opt.neqn_p,1);

if strcmp(study.solver,'Direct')==1

    tic;
    spparms('spumoni',1)
    U = opt.LHS \ (opt.RHS);
    time = toc;
    fprintf('Solution found in %f s\n',time);

    opt.u1 = U(1:opt.neqn_u);
    opt.u2 = U(opt.neqn_u+1:2*opt.neqn_u);
    opt.p = U(2*opt.neqn_u+1:end);

elseif strcmp(study.solver,'pcg')==1
    
    tol = 1e-6;
    opt.p = pcg(S, RHS, [], [], opt.M);
    opt.u1 = A1 \ (D1' * opt.p + Bf1);
    opt.u2 = A2 \ (D2' * opt.p + Bf2);

elseif strcmp(study.solver,'Uzawa')==1
    max_iter = 50;
    tol = 1e-12;
    [opt.p, conv_message] = pcgUzawa(S, Bf1, Bf2, Fp, opt.M, p0, A1, A2, D1, D2, tol, max_iter);
    opt.u1 = A1 \ (D1' * opt.p + Bf1);
    opt.u2 = A2 \ (D2' * opt.p + Bf2);
disp(conv_message)
elseif strcmp(study.solver,'UzawamodJGA')==1

    [opt.p] = pcg_mod(S, Bf1, Bf2, Fp, opt.M, p0, A1, A2, D1, D2);
    opt.u1 = A1 \ (D1' * opt.p + Bf1);
    opt.u2 = A2 \ (D2' * opt.p + Bf2);

end

sum(abs(A1*opt.u1-D1'*opt.p-Bf1))
sum(abs(A2*opt.u2-D2'*opt.p-Bf2))
sum(abs(-D1*opt.u1-D2*opt.u2-Fp))

end
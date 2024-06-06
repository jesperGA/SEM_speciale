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

H1 = opt.LHS_BC(1:opt.neqn_u,1:opt.neqn_u);
H2 = opt.LHS_BC(opt.neqn_u+1:2*opt.neqn_u,opt.neqn_u+1:2*opt.neqn_u);
D1 = -opt.LHS_BC(2*opt.neqn_u+1:end,1:opt.neqn_u);
D2 = -opt.LHS_BC(2*opt.neqn_u+1:end,opt.neqn_u+1:2*opt.neqn_u);

S = D1 * (H1 \ D1') + D2 * (H2 \ D2');

dt = study.t(2) - study.t(1);

mu = mesh.Material(1);
rho = mesh.Material(3);

opt.u1 = sparse(opt.neqn_u,length(study.t));
opt.u2 = sparse(opt.neqn_u,length(study.t));
opt.p  = sparse(opt.neqn_p,length(study.t));

u0=sparse(opt.neqn,1);
u0(find(opt.g)) =opt.g(find(opt.g));
opt.u1(:,1) = u0(1:opt.neqn_u);
opt.u2(:,1) = u0(opt.neqn_u+1:2*opt.neqn_u);
opt.p(:,1) = u0(2*opt.neqn_u+1:end);

MAXIT = 1e2;
tol = 1e-6;

E = opt.D1 * (opt.B \ opt.D1') + opt.D2 * (opt.B \ opt.D2'); %OBS w/o BCs enforced
preCond = (mu * inv(opt.M) + rho/dt * inv(E)); % Lav evt. noget chol() el. lign.


for n = 2:length(study.t)

    g1 = opt.B * (opt.f_u1 + rho / dt * opt.u1(:,n-1));
    g2 = opt.B * (opt.f_u2 + rho / dt * opt.u2(:,n-1));

    RHS = [g1; g2; sparse(opt.neqn_p,1)];

    [~, RHS] = applyBoundaryConditions(opt.LHS, RHS, opt.Null, opt.g, opt.neqn);
    RHS1 = RHS(1:opt.neqn_u);
    RHS2 = RHS(opt.neqn_u+1:2*opt.neqn_u);
    RHS3 = RHS(2*opt.neqn_u+1:end);
    RHS_Uzawa = - D1 * (H1 \ RHS1) - D2 * (H2 \ RHS2) - RHS3;

    if strcmp(study.solver,'Direct')

        tic;
        U = opt.LHS_BC \ (RHS);
        time = toc;
        fprintf('Solution found in %f s\n',time);

        opt.u1(:,n) = U(1:opt.neqn_u);
        opt.u2(:,n) = U(opt.neqn_u+1:2*opt.neqn_u);
        opt.p(:,n)= U(2*opt.neqn_u+1:end);
        opt.p(:,n) = opt.p(:,n)-mean(opt.p(:,n));

    else
        if strcmp(study.solver,'pcg')
            if n==2;
                preCond = inv(preCond);
            end
            opt.p(:,n) = pcg(S, RHS_Uzawa, tol, MAXIT, preCond, [] , opt.p(:,n-1));
    
        elseif strcmp(study.solver,'Uzawa')
    
            opt.p(:,n) = pcgUzawa(S, RHS_Uzawa, preCond,  opt.p(:,n-1), H1, H2, D1, D2, tol, MAXIT);
    
        elseif strcmp(study.solver,'UzawamodJGA')
    
            opt.p(:,n) = pcg_mod(S, RHS_Uzawa, preCond, opt.p(:,n-1), H1, H2, D1, D2, tol, MAXIT);
    
        end
        opt.u1(:,n) = H1 \ (D1' * opt.p(:,n) + RHS1);
        opt.u2(:,n) = H2 \ (D2' * opt.p(:,n) + RHS2);
    end


end

sum(abs(H1*opt.u1(:,n)-D1'*opt.p(:,n)-RHS1))
sum(abs(H2*opt.u2(:,n)-D2'*opt.p(:,n)-RHS2))
sum(abs(-D1*opt.u1(:,n)-D2*opt.u2(:,n)-RHS3))

end
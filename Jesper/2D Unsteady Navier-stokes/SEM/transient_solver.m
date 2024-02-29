function [opt] = transient_solver(opt,study,mesh)

pref_dof = mesh.pref_dof;
dt = study.dt;
ndt = study.nt;

xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);

if strcmp(study.p_type,'roenquist') == 1
    opt.P1 = 2+pi*cos(pi*xx).*sin(pi*yy);
    opt.P2 = pi*sin(pi*xx).*cos(pi*yy);
elseif strcmp(study.p_type,'bercover') == 1

    F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
        2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2))+X2-1/2;
    F2 = @(X1,X2)   F1(X2,X1);

    opt.P1 = F1(xx,yy);
    opt.P2 = F2(xx,yy);
else
    disp('No load defined')
    opt.P1 = zeros(length(xx),1);
    opt.P2 = zeros(length(yy),1);
end

if strcmp(study.int_type,'BDFk') == 1 || strcmp(study.int_type,'BDF1AB3') == 1
    gamma = [1, 1 , 1];
    b0 = gamma(2)/gamma(1);
    a_k = gamma(3:end)./gamma(1);
    %Only for k=1 order in BDF-k [b0 = gamma(2)/gamma(1), a_k =gamma(3:end)./gamma(1)];


    % DEFINE BETA VARIABLES FOR BDF.
    beta = [flip(a_k)/b0, 1/b0];

    AB_facs = [23/12, -4/3, 5/12];

    H = beta(end)/dt*opt.M+opt.K;

    %SET UP GLOBAL MATRICES FOR BCs.
    null_sys = [opt.Null1, zeros(opt.neqnV),zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnV),opt.Null2,zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),speye(opt.neqnP,opt.neqnP)];
    g_sys =  [opt.g1;opt.g2;zeros(opt.neqnP,1)];


    %Assembly LHS(sys_mat) and part of RHS (B_mat) of the direct equation
    sys_mat = [H,zeros(opt.neqnV,opt.neqnV),-opt.DE1.';
        zeros(opt.neqnV,opt.neqnV), H,-opt.DE2.';
        -opt.DE1,-opt.DE2,zeros(opt.neqnP,opt.neqnP)];
    B_mat = [opt.M,zeros(opt.neqnV,opt.neqnV),zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnV,opt.neqnV),opt.M,zeros(opt.neqnV,opt.neqnP);
        zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnP)];


    free = diag(null_sys);
    sys_org = sys_mat;

    sys_mat = null_sys'*sys_mat*null_sys-(null_sys-speye(size(null_sys)));
    opt.sys_mat = sys_mat;

    BCs = sys_org*g_sys;

    U(:,1) = study.U0;
    opt.Pr(:,1) = zeros(opt.neqnP,1);

    [L,Up,Pp] = lu(sys_mat);

    for i = 2:ndt
        Un = [];
        Un = [U;zeros(opt.neqnP,1)];
        if strcmp(study.int_type,'BDFk') == 1
            P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]+(beta(end-1)/dt)*Un);

        else
            J = min(i,3); %Define order of Adam bashforth. Will make sure we dont exceed array limits.
            CV = 0;
            for j=1:J
                C_temp = adv_mat_assembly(mesh,opt,study,U);
                CV = CV+AB_facs(j)*C_temp;

            end

            P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]+(beta(end-1)/dt)*Un)-study.RE*CV;

        end

        P = P-BCs;
        P(~free) = 0;
        P(find(g_sys)) = g_sys(find(g_sys));

        if strcmp(study.solve_type,'direct') == 1
            if strcmp(study.direct_type,'LU') == 1
                y = L\(Pp*P);
                sol = Up\y;
            elseif strcmp(study.direct_type,'backslash') == 1
                sol = sys_mat \ (P);
            else
                y = L\(Pp*P);
                sol = Up\y;
            end
            opt.Pr(:,i) = sol(end-opt.neqnP+1:end);
            opt.Pr(:,i) = opt.Pr(:,i)-mean(opt.Pr(:,i));
            opt.U(:,i) = sol(1:2*opt.neqnV);
            U = sol(1:2*opt.neqnV);

        elseif strcmp(study.solve_type,'uzawa') == 1

            neqnp = opt.neqnP;
            neqnv = opt.neqnV;

            D1 = -sys_mat(1:neqnv,2*neqnv+1:end)';
            D2 = -sys_mat(neqnv+1:2*neqnv,2*neqnv+1:end)';
            H1 = sys_mat(1:neqnv,1:neqnv);
            H2 = sys_mat(neqnv+1:2*neqnv,neqnv+1:2*neqnv);

            RHS1 = P(1:neqnv);RHS2 = P(neqnv+1:neqnv*2);RHS3 = P(neqnv*2+1:end);
            A = D1*(H1\D1')+D2*(H2\D2');
            B = -D1*(H1\RHS1)-D2*(H2\RHS2)-RHS3;

            A = (A+A')/2;

            if strcmp(study.precon,'mhat') == 1

                p = pcg_mod(A,B,opt.Mh,opt.Pr(:,i-1),H1,H2,D1,D2);

            elseif strcmp(study.precon,'P') == 1

                E = opt.DE1*inv(opt.M)*opt.DE1'+opt.DE2*inv(opt.M)*opt.DE2';
                Pc = inv(opt.Mh) + 1/dt*inv(E);

                p = pcg_mod(A,B,inv(Pc),opt.Pr(:,i-1),H1,H2,D1,D2);
                % p = pcg(A,B,1e-6,100,Pc,[],opt.Pr(:,i-1));
            end

            opt.Pr(:,i) = p;
            U = [H1 \ (D1'*p + RHS1)
                H2 \ (D2'*p + RHS2)];
            opt.U(:,i) = U;


        end


    end

else
    disp('Time integration scheme not implemented')
end
end
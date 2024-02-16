function [opt] = solver(opt, study,mesh)
xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
xxp = mesh.Xp(:,2);yyp = mesh.Xp(:,3);
if strcmp(study.p_type,'roenquist') == 1
    opt.P1 = 2+pi*cos(pi*xx).*sin(pi*yy);
    opt.P2 = pi*sin(pi*xx).*cos(pi*yy);
elseif strcmp(study.p_type,'bercover') == 1
    opt.P1 = (1536 .* xx .^ 2 .* (xx - 1) .^ 2 .* (2 .* yy - 1)) + (256 .*...
        (yy - 1) .* (2 .* yy - 1) .* yy .* (12 .* xx .^ 2 - 12 .* xx + 2));
    opt.P2 = (1536 .* yy .^ 2 .* (yy - 1) .^ 2 .* (2 .* xx - 1)) + (256 .*...
        (xx - 1) .* (2 .* xx - 1) .* xx .* (12 .* yy .^ 2 - 12 .* yy + 2));
else
    error('Load case not defined')
end

%ENFRORCE BOUNRDARIES enforce somehow else?
% if strcmp(study.bound_type,'total') == 1
%Setup global matrix and vector for BCS
null_sys = [opt.Null1, zeros(opt.neqnV),zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnV),opt.Null2,zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),speye(opt.neqnP,opt.neqnP)];
g_sys =  [opt.g1;opt.g2;zeros(opt.neqnP,1)];


%Assembly LHS(sys_mat) and RHS (B_mat*opt.P) of the direct equation
sys_mat = [opt.K,zeros(opt.neqnV,opt.neqnV),-opt.DE1.';
    zeros(opt.neqnV,opt.neqnV), opt.K,-opt.DE2.';
    -opt.DE1,-opt.DE2,zeros(opt.neqnP,opt.neqnP)];
B_mat = [opt.M,zeros(opt.neqnV,opt.neqnV),zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnV,opt.neqnV),opt.M,zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnP)];
opt.P = B_mat*[opt.P1;opt.P2;zeros(opt.neqnP,1)];

free = diag(null_sys);
sys_org = sys_mat;
opt.P = opt.P-sys_org*g_sys;
% opt.P
opt.P(~free) = 0;
opt.P(find(g_sys)) = g_sys(find(g_sys));
% opt.P
sys_mat = null_sys'*sys_mat*null_sys-(null_sys-speye(size(null_sys)));
% elseif strcmp(study.bound_type,'local') == 1
opt.sys_mat = sys_mat;

% end

if strcmp(study.solve_type,'direct') == 1

    opt.U = sys_mat \ (opt.P);

elseif strcmp(study.solve_type,'uzawa') == 1
    neqnp = opt.neqnP;
    neqnv = opt.neqnV;

    D1 = -sys_mat(1:neqnv,2*neqnv+1:end)';
    D2 = -sys_mat(neqnv+1:2*neqnv,2*neqnv+1:end)';
    K1 = sys_mat(1:neqnv,1:neqnv);
    K2 = sys_mat(neqnv+1:2*neqnv,neqnv+1:2*neqnv);

    RHS1 = opt.P(1:neqnv);RHS2 = opt.P(neqnv+1:neqnv*2);RHS3 = opt.P(neqnv*2+1:end);
    A = D1*(K1\D1')+D2*(K2\D2');
    B = -D1*(K1\RHS1)-D2*(K2\RHS2)-RHS3;

    A = (A+A')/2;

    p = pcg_mod(A,B,opt.Mh,zeros(opt.neqnP,1),K1,K2,D1,D2);
    % [p,pcg_flag] = pcg(A,B,[],[],opt.Mh);
    % if  pcg_flag ~=0
    %     warning(['Uzawa: PCG not converged succesfully. Flag: ',num2str(pcg_flag)])
    % end

    opt.U = [K1\D1'*p+K1\RHS1
        K2\D2'*p+K2\RHS2
        p];
end

end
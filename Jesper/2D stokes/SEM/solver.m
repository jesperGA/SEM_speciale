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
end
if strcmp(study.solve_type,'direct') == 1
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
    opt.P(~free) = full(g_sys(~free));
    % opt.P
    sys_mat = null_sys'*sys_mat*null_sys-(null_sys-speye(size(null_sys)));
    % B_mat = null_sys'*B_mat*null_sys-(null_sys-speye(size(null_sys)));
    % % Solve static problem
    % (B_mat*opt.P)
    opt.U = sys_mat \ (opt.P);

elseif strcmp(study.solve_type,'uzawa') == 1
    
end

end
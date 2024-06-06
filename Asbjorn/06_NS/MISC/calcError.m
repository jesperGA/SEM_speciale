function [error_u, error_p] = calcError(mesh,study,opt)
    x = mesh.X(:,2);
    y = mesh.X(:,3);
    xp = mesh.Xp(:,2);
    yp = mesh.Xp(:,3);
    U1 = study.U1;
    U2 = study.U2;
    P = study.P;

    error_u_t = zeros(length(study.t),1);
    error_p_t = zeros(length(study.t),1);
    for n = length(study.t_steps):length(study.t_steps)
        time = study.t(study.t_steps(n));
        error_u_t(n) = norm([opt.u1(:,n);opt.u2(:,n)] - [U1(x,y,time);U2(x,y,time)],1);
        error_p_t(n) = norm(opt.p(:,n) - P(xp,yp,time),1);
    end
    error_u = norm(error_u_t,'inf');
    error_p = norm(error_p_t,'inf');
end
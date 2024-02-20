function [error_u, error_p] = calcError(mesh,study,opt)
    x = mesh.X(:,2);
    y = mesh.X(:,3);
    U1 = study.U1;
    U2 = study.U2;
    P = study.P;
    error_u = norm([opt.u1;opt.u2] - [U1(x,y);U2(x,y)],'inf');
    x = mesh.Xp(:,2);
    y = mesh.Xp(:,3);
    error_p = norm(opt.p - P(x,y),'inf');
end
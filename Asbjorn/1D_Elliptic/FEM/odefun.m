function [z_p] = odefun(t, z, M, C, K, tvec, fvec)
    
    neqn=size(M,1);
    z_p = zeros(2*neqn,1);
    z_p(1:neqn) = z(neqn+1:end);
    fval = interp1(tvec, fvec, t);
    f=zeros(neqn,1);
    f(1)=fval;
    z_p(neqn+1:2*neqn) = M \ (f - C*z(neqn+1:end) -K*z(1:neqn));
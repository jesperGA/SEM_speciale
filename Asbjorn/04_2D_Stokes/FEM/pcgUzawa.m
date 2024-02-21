function [p] = pcgUzawa(S, RHS, M, p0, A1, A2, D1, D2)
max_iter = 20;
tol = 1e-12;

u = p0;

r_old = - RHS + S*p0;
q_old= M \ r_old;
p = q_old;

for m = 1:50

    y1 = D1'*p;
    y2 = D2'*p;

    z1 = A1\y1; z2 = A2\y2;

    Sp = D1*z1+D2*z2;

    a = -(dot(q_old,r_old)/dot(p,Sp));
    u = u + a*p;
    r = r_old + a*Sp;
    q = M\r;
    b = dot(q,r)/dot(q_old,r_old);
    p = q + b * p; % <- wtf, fejl i Rønquist?

    ch = norm(r);

    r_old = r; q_old = q;

    if i>=max_iter || ch < tol
        p = u;
        if i>=max_iter

            warning('pcg not converged succesfully. Max iter reached')
        end

        break
    end
end

end
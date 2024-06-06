function [p,conv_message] = pcgUzawa(S, RHS, M, p0, A1, A2, D1, D2, tol, max_iter)

u = p0;

r_old = - RHS + S*p0;
q_old= M * r_old;
p = q_old;

for m = 1:max_iter

    y1 = D1'*p;
    y2 = D2'*p;

    z1 = A1\y1; z2 = A2\y2;

    Sp = D1*z1+D2*z2;

    a = -(dot(q_old,r_old)/dot(p,Sp));
    u = u + a*p;
    r = r_old + a*Sp;
    q = M*r;
    b = dot(q,r)/dot(q_old,r_old);
    p = q + b * p; % <- wtf, fejl i Rønquist?

    ch = norm(r);

    r_old = r; q_old = q;

    if m>=max_iter || ch < tol
        p = u;
        if m>=max_iter

            warning('pcg not converged succesfully. Max iter reached')
        end

        break
    end
end

conv_message = itermsg('pcg',tol,max_iter,m,0,m,ch);

end
function [p] = pcg_mod(A,B,Mh,p0,K1,K2,DE1,DE2)
max_iter = 50;
tol = 1e-6;

x_old = p0;

r_old = B-A*p0;
z_old = Mh\r_old;
w_old = z_old;

for i =1:max_iter

    y1 = DE1'*w_old;
    y2 = DE2'*w_old;

    z1 = K1\y1; z2 = K2\y2;

    sub = DE1*z1+DE2*z2;

    alpha = dot(r_old,z_old)/(dot(sub,w_old));
    x = x_old+alpha*w_old;
    r = r_old-alpha*(sub);
    z = Mh\r;
    beta = dot(r,z)/dot(r_old,z_old);
    w = z+beta*w_old;

    ch = norm(r);
    
    x_old = x; r_old = r; z_old = z; w_old = w;
    
    if i>=max_iter || ch < tol
        p = x;
        if i>=max_iter
            
            warning('pcg not converged succesfully. Max iter reached')
        end

        break
    end
end

end
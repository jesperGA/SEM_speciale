function [Edis] = energyDissipated(mu,u1,u2,grad1,grad2,Alpha,w,J)

    du1dx1 = grad1 * u1;
    du2dx2 = grad2 * u2;
    du2dx1 = grad1 * u2;
    du1dx2 = grad2 * u1;

    Edis_nodal = mu * (2 * du1dx1.^2 + 2 * du2dx2.^2 + (du2dx1 + du1dx2).^2) + Alpha * (u1.^2 + u2.^2);

    Edis_nodal = reshape(Edis_nodal,length(w),length(w));

    Edis = 0;
    for i = 1:length(w)
        for j = 1:length(w)
            Edis = Edis + Edis_nodal(i,j) * w(i) * w(j) * J(i,j);
        end
    end

end

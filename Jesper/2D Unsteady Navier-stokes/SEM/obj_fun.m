function [phi,dissE] = obj_fun(opt,T,mu,IX,X,w)

D1 = opt.D1;
D2 = opt.D2;




if isfield(opt,'alpha')
    alpha = opt.alpha;
else 
    alpha = zeros(size(D1));
end
%% Global matrices
% outcomment for global matrices
% U1 = opt.U(1:opt.neqnV,:);
% U2 = opt.U(opt.neqnV+1:end,:);
% D = [D1 zeros(size(D2))
%       zeros(size(D2)) D2];
% U = [U1;U2];

% A = [alpha zeros(size(D1))
%     zeros(size(D1)) alpha];
%%
diss_elemental = zeros([opt.nel,numel(T)]);

for i = 1:numel(T)
    for e = 1:opt.nel
        
        nen = IX(:,:,e);
        xy = X(nen,2:3);
        edof = nen(:);

        u1 = U1(edof);
        u2 = U2(edof);
        d1 = D1(edof,edof);
        d2 = D2(edof,edof);
        alpha_loc = alpha(edof,edof);

        [diss_elemental(e,i)] = energyDissipated(mu,u1,u2,d1,d2,alpha_loc,w,opt.dJ(:,:,e));
    end
    dissE = sum(diss_elemental,1);
    % dissE(i) = 1/2*(mu*((D*U(:,i)).'*(D*U(:,i))+(A*U(:,i)).'*(A*U(:,i))));
end

phi = trapz(T,dissE);

end
%% Energy calculator per element
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
% function [Cv] = adv_mat_assembly2(IXv, ME_DE1, ME_DE2,n_GLL,neqn,U,S)
function [Cv] = adv_mat_assembly2(IXv, ME_DE1, ME_DE2, n_GLL, neqn, U)
neqnV = neqn(1);
neqnP = neqn(2);

nel = size(ME_DE1, 3);
ldof = n_GLL^2;

U1 = U(1:neqnV);
U2 = U(neqnV+1:end);

resultV1 = zeros(neqnV, 1);
resultV2 = zeros(neqnV, 1);

% localContributionV1{nel}=[];
% localContributionV2{nel}=[];

for e = 1:nel
    nen = IXv(:,:,e);
    nen = nen(:);

    v1_vec = U1(nen);
    v2_vec = U2(nen);

    v1 = sparse(1:ldof,1:ldof,v1_vec,ldof,ldof);
    v2 = sparse(1:ldof,1:ldof,v2_vec,ldof,ldof);



    ME1 = ME_DE1(:,:,e);
    ME2 = ME_DE2(:,:,e);

    % Compute the local contributions directly
    ce = v1*ME1 +v2* ME2;
    localContributionV1 = ce*v1_vec;
    localContributionV2 = ce*v2_vec;
    % localContributionV1 = (ME1 * v1 + ME2 * v2)*v1_vec;
    % localContributionV2 = (ME1 * v1 + ME2 * v2)*v2_vec;
    % accumulated_indices = sub2ind([neqnV, 1], nen, nen);

    % Accumulate the contributions into the global result vectors
    resultV1(nen) = resultV1(nen) + localContributionV1;
    resultV2(nen) = resultV2(nen) + localContributionV2;
end
% for e = 1:nel
%     nen = IXv(:,:,e);
%     nen = nen(:);
%     resultV1(nen) = resultV1(nen) + localContributionV1{e};
%     resultV2(nen) = resultV2(nen) + localContributionV2{e};
% end

% Combine results into the final vector
Cv = [resultV1;
    resultV2;
    sparse(neqnP, 1)];  % Presuming that no pressure coupling is computed here
end
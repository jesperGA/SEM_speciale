function [opt] = transient_solver(opt,study,mesh)

% pref_dof = mesh.pref_dof;
dt = study.dt;
ndt = study.nt;

xx = mesh.Xv(:,2);yy = mesh.Xv(:,3);
neqnV = opt.neqnV;

if strcmp(study.p_type,'roenquist') == 1
    opt.P1 = 2+pi*cos(pi*xx).*sin(pi*yy);
    opt.P2 = pi*sin(pi*xx).*cos(pi*yy);
elseif strcmp(study.p_type,'bercover') == 1

    F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
        2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2))+X2-1/2;
    F2 = @(X1,X2)   F1(X2,X1);

    opt.P1 = F1(xx,yy);
    opt.P2 = F2(xx,yy);
else
    disp('No load defined')
    opt.P1 = zeros(length(xx),1);
    opt.P2 = zeros(length(yy),1);
end


gamma = [1, 1 , 1];
b0 = gamma(2)/gamma(1);
a_k = gamma(3:end)./gamma(1);
%Only for k=1 order in BDF-k [b0 = gamma(2)/gamma(1), a_k =gamma(3:end)./gamma(1)];

%Define integration variables
if strcmp(study.int_type,'BDFk') == 1 || strcmp(study.int_type,'BDF1AB3') == 1 || strcmp(study.int_type,'explicit') == 1
    % DEFINE BETA VARIABLES FOR BDF.
    beta = [flip(a_k)/b0, 1/b0];
    AB_facs = [23/12, -4/3, 5/12];
    if strcmp(study.int_type,'explicit') == 1
        if isfield(mesh,'material') && length(mesh.material) >= 2
            rho = mesh.material(1); mu = mesh.material(2);
            H = rho*beta(end)/dt*opt.M;
        else
            H = beta(end)/dt*opt.M;
        end
    else
        if isfield(mesh,'material') && length(mesh.material) >= 2
            rho = mesh.material(1); mu = mesh.material(2);
            H = rho*beta(end)/dt*opt.M+mu*opt.K;
        else
            H = beta(end)/dt*opt.M+opt.K;
        end
    end
elseif strcmp(study.int_type,'BDF3EX3') == 1
    BF3_0 = 11/6;
    BF3_fac = [3,-3/2, 1/3];
    EX3_fac = [3,-3,1];
    if isfield(mesh,'material') && length(mesh.material) >= 2
        rho = mesh.material(1); mu = mesh.material(2);
        H = rho*BF3_0*(1/dt)*opt.M+mu*opt.K;
    else
        H = BF3_0*(1/dt)*opt.M+opt.K;
    end
end

%Add brinkmann term
if isfield(study,'Brinkmann') && study.Brinkmann
    body_ind = mesh.body_ind(:);
    alpha_vec = study.alpha;
    alpha = sparse(body_ind,body_ind,mesh.alpha_fac(:,2)*alpha_vec,opt.neqnV,opt.neqnV);

    H = H+alpha;
    opt.alpha = alpha;
end

% Precompute matrices for C.
ldof = (study.n_GLL)^2;
ME_DE1 = zeros(size(opt.DE1v));
ME_DE2 = zeros(size(opt.DE2v));

I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);

ntriplets = 0;
for e =1:opt.nel
    nen = mesh.IXv(:,:,e);
    nen = nen(:);

    edof = zeros(ldof,1);
    indices = 1:length(nen);
    edof(indices) = nen;

    % Create all combinations of krow and kcol
    % [krow, kcol] = meshgrid(1:ldof, 1:ldof);
    % krow = krow(:);
    % kcol = kcol(:);

    % Vectorize the assignment to I, J, CE
    % I(ntriplets+1:ntriplets+ldof^2) = edof(krow);
    % J(ntriplets+1:ntriplets+ldof^2) = edof(kcol);

    ME_DE1(:,:,e) = opt.ME(:,:,e) * opt.DE1v(:,:,e);
    ME_DE2(:,:,e) = opt.ME(:,:,e) * opt.DE2v(:,:,e);

    ntriplets = ntriplets + ldof^2;

end
% S_template = sparse(I, J, zeros(size(I)), max(I), max(J));
% [I, J, ME_DE1_flat, ME_DE2_flat] = precomputeGlobalCmatr(ldof, opt.nel,opt.ME , opt.DE1v,opt.DE2v, mesh.IXv);
% [I, J, ME_DE1_flat, ME_DE2_flat] = precomputeCComponents(study, mesh, opt);
%%
%SET UP GLOBAL MATRICES FOR BCs.
null_sys = [opt.Null1, zeros(opt.neqnV),zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnV),opt.Null2,zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),speye(opt.neqnP,opt.neqnP)];


%Assembly LHS(sys_mat) and part of RHS (B_mat) of the direct equation
sys_mat = [H,zeros(opt.neqnV,opt.neqnV),-opt.DE1.';
    zeros(opt.neqnV,opt.neqnV), H,-opt.DE2.';
    -opt.DE1,-opt.DE2,zeros(opt.neqnP,opt.neqnP)];
B_mat = [opt.M,zeros(opt.neqnV,opt.neqnV),zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnV,opt.neqnV),opt.M,zeros(opt.neqnV,opt.neqnP);
    zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnV),zeros(opt.neqnP,opt.neqnP)];


free = diag(null_sys);
sys_org = sys_mat;

sys_mat = null_sys'*sys_mat*null_sys-(null_sys-speye(size(null_sys)));
opt.sys_mat = sys_mat;
if strcmp(study.BC_type,'static') == 1
    g_sys =  [opt.g1;opt.g2;zeros(opt.neqnP,1)];
    BCs = sys_org*g_sys;
else
    g_sys = [opt.g_sys(xx,yy,0);zeros(opt.neqnP,1)];
    g_sys(free==1) = 0;
end

if strcmp(study.int_type,'explicit') == 1

    D1 = -sys_mat(1:neqnV,2*neqnV+1:end)';
    D2 = -sys_mat(neqnV+1:2*neqnV,2*neqnV+1:end)';
    H1 = sys_mat(1:neqnV,1:neqnV);
    H2 = sys_mat(neqnV+1:2*neqnV,neqnV+1:2*neqnV);

    sys_mat =-1*(D1*(H1\D1')+D2*(H2\D2'));

end
frames = ceil(study.fps*study.T)+1;
fdt = 1/study.fps;

opt.U = zeros(opt.neqnV*2,frames);
opt.P = zeros(opt.neqnP,frames);

%%
fname = study.name;
folderName = [fname, '_anis'];  % Construct folder name
if ~exist(folderName, 'dir')
    mkdir(folderName);          % Create the folder if it doesn't exist
end
file_path = fullfile(folderName, [fname, '.nek5000']);

% Open the file for writing with .nek5000 extension
fileID = fopen(file_path, 'w+');

% Write the specified lines with the variable contents into the file
fprintf(fileID, 'filetemplate: %s%%01d.f%%05d\n', fname);
fprintf(fileID, 'firsttimestep: %d\n', 1);
fprintf(fileID, 'numtimesteps: %d\n', frames);

% Close the file
fclose(fileID);

%%

U = study.U0;
opt.Pr(:,1) = zeros(opt.neqnP,1);
opt.U(:,1) = U;
disp('Factorizing system matrix')
tic;
% [L,Up,Pp] = lu(sys_mat);
dsys_mat = decomposition(sys_mat,'banded');

% dsys_mat = decomposition(sys_mat);
fprintf('Done decomposotioning. Time: %2.2f seconds\n Decomp type = %s\n',toc,dsys_mat.Type);

%Prepare for CFL calculations:
[xLength, yLength] = findNearestXYDistances(mesh.Xv(:,2:3));
closestDeltas = [xLength,yLength];
pre_cfl = dt./closestDeltas;


% fprintf('Solver status:') %For AeStheTics
prev_mes = 0;
prev_time = 0;
fprintf('Solver status: ')

%For advection matrix
IXv = mesh.IXv;ME = opt.ME;DE1v = opt.DE1v;DE2v = opt.DE2v; n_GLL = study.n_GLL;
Cold = sparse(2*opt.neqnV+opt.neqnP,1);
Coldold  = Cold;

Unold = zeros(2*opt.neqnV+opt.neqnP,1);Unoldold = Unold;
Unold(~free) = g_sys(~free);Unoldold(~free) = g_sys(~free);

count = 1;
opt.movie_step(count) = 1;
save_step(folderName,fname,mesh,study.N_interp,opt.U(1:2*opt.neqnV)...
    ,opt.Pr(:,count),opt.neqnV,study.t(1),1,count,opt.nel)
next_frame_time = fdt;
for i = 2:ndt
    %% For AeStheTics
    status_procent = (i/ndt)*100;
    if status_procent > prev_time+0.5
        mess = sprintf('%3.2f percent',status_procent);
        back_space = repmat('\b',1,prev_mes);
        fprintf(back_space);
        fprintf(mess);
        prev_mes = length(mess);
        prev_time = status_procent;

    end
    %%
    Un = [U;zeros(opt.neqnP,1)];
    if strcmp(study.int_type,'BDFk') == 1
        P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]+(beta(end-1)/dt)*Un);

    elseif strcmp(study.int_type,'BDF1AB3') == 1 || strcmp(study.int_type,'explicit') == 1
        % C = adv_mat_assembly_mex(IXv,ME,DE1v,DE2v,n_GLL,[opt.neqnV,opt.neqnP],opt.U(:,i-1));
        % C_= adv_mat_assembly(IXv,ME,DE1v,DE2v,n_GLL,[opt.neqnV,opt.neqnP],opt.U(:,i-j));
        C = adv_mat_assembly2(IXv, ME_DE1, ME_DE2,n_GLL,[opt.neqnV,opt.neqnP],U(1:2*opt.neqnV));
        % C3 = adv_mat_assembly(IXv, ME_DE1, ME_DE2,n_GLL,[opt.neqnV,opt.neqnP],Un(1:2*opt.neqnV),S_template);
        % C = assembleC(I, J, ME_DE1_flat, ME_DE2_flat, U, opt.neqnV, opt.neqnP);

        CV = AB_facs(1)*C+AB_facs(2)*Cold+AB_facs(3)*Coldold;
        Coldold = Cold; Cold = C;
        if strcmp(study.int_type,'BDF1AB3') == 1
            if isfield(mesh,'material') && length(mesh.material)>=2
                P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]+(rho*beta(end-1)/dt)*Un)-rho*CV;
            else
                P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]+(beta(end-1)/dt)*Un)-study.RE*CV;
            end
        elseif strcmp(study.int_type,'explicit') == 1
            u1 = U(1:neqnV); u2 = U(neqnV+1:end);
            P1 = -mu * opt.K * u1 + opt.M * (rho / dt * (beta(end-1) * u1));
            P2 = -mu * opt.K * u2 + opt.M * (rho / dt * (beta(end-1) * u2));
            P = [P1;P2;zeros(opt.neqnP,1)] - rho*CV;

        end

    elseif strcmp(study.int_type,'BDF3EX3') == 1
        C = adv_mat_assembly2(IXv, ME_DE1, ME_DE2,n_GLL,[opt.neqnV,opt.neqnP],U);
        % C2 = adv_mat_assembly_mex(IXv,ME,DE1v,DE2v,n_GLL,[opt.neqnV,opt.neqnP],Un(1:2*opt.neqnV));
        % C_= adv_mat_assembly(IXv,ME,DE1v,DE2v,n_GLL,[opt.neqnV,opt.neqnP],opt.U(:,i-j));
        % C = globalCmatr_ADJ(I, J, ME_DE1, ME_DE2, opt.nel, U(1:opt.neqnV), U(opt.neqnV+1:end), ldof, [opt.neqnV,opt.neqnP]);
        % C = assembleC(I, J, ME_DE1_flat, ME_DE2_flat, U, opt.neqnV, opt.neqnP);
        CV = EX3_fac(1)*C+EX3_fac(2)*Cold+EX3_fac(3)*Coldold;
        Coldold = Cold; Cold = C;


        if isfield(mesh,'material') && length(mesh.material)>=2
            P = B_mat*([opt.P1;opt.P2;zeros(opt.neqnP,1)]+(rho/dt)*(BF3_fac(1) * Un + BF3_fac(2) * Unold + BF3_fac(3) * Unoldold))-rho*CV;
        else
            P = B_mat * ( [opt.P1;opt.P2;zeros(opt.neqnP,1)] + (1/dt) * (BF3_fac(1) * Un + BF3_fac(2) * Unold + BF3_fac(3) * Unoldold)) - study.RE*CV;
        end
        Unoldold = Unold;Unold = Un;
    end


    if strcmp(study.BC_type,'dynamic') == 1
        g_sys = [opt.g_sys(xx,yy,study.t(i));zeros(opt.neqnP,1)];
        g_sys(free==1) = 0;
        BCs = sys_org*g_sys;
    end
    %Prepare RHS with BCs
    P = P-BCs;
    P(~free) = 0;
    P(find(g_sys)) = g_sys(find(g_sys));
    if strcmp(study.int_type,'explicit') == 1
        RHS1 = P(1:opt.neqnV);RHS2 = P(opt.neqnV+1:2*opt.neqnV);
        RHS3 = P(2*opt.neqnV+1:end);
        p = dsys_mat\(D1*(H1\RHS1)+D2*(H2\RHS2)+RHS3);
        U = [H1\(RHS1+D1'*p);
            H2\(RHS2+D2'*p)];
        sol = [U;p];

        CFL = U(1:neqnV).*pre_cfl(:,1)+U(neqnV+1:end).*pre_cfl(:,2);
        opt.CFL_m(i) = max(CFL);
    else
        if strcmp(study.solve_type,'direct') == 1
            if strcmp(study.direct_type,'LU') == 1
                y = L\(Pp*P);
                sol = Up\y;
            elseif strcmp(study.direct_type,'backslash') == 1
                sol = dsys_mat \ (P);
            else
                y = L\(Pp*P);
                sol = Up\y;
            end



            U = sol(1:2*opt.neqnV);

        elseif strcmp(study.solve_type,'uzawa') == 1

            neqnp = opt.neqnP;
            neqnv = opt.neqnV;

            D1 = -sys_mat(1:neqnv,2*neqnv+1:end)';
            D2 = -sys_mat(neqnv+1:2*neqnv,2*neqnv+1:end)';
            H1 = sys_mat(1:neqnv,1:neqnv);
            H2 = sys_mat(neqnv+1:2*neqnv,neqnv+1:2*neqnv);

            RHS1 = P(1:neqnv);RHS2 = P(neqnv+1:neqnv*2);RHS3 = P(neqnv*2+1:end);
            A = D1*(H1\D1')+D2*(H2\D2');
            B = -D1*(H1\RHS1)-D2*(H2\RHS2)-RHS3;

            A = (A+A')/2; %Enforce strict symmetry

            if strcmp(study.precon,'mhat') == 1

                p = pcg_mod(A,B,opt.Mh,opt.Pr(:,i-1),H1,H2,D1,D2);

            elseif strcmp(study.precon,'P') == 1

                E = opt.DE1*inv(opt.M)*opt.DE1'+opt.DE2*inv(opt.M)*opt.DE2';
                Pc = inv(opt.Mh) + 1/dt*inv(E);

                p = pcg_mod(A,B,inv(Pc),opt.Pr(:,i-1),H1,H2,D1,D2);
                % p = pcg(A,B,1e-6,100,Pc,[],opt.Pr(:,i-1));
            end

            opt.Pr(:,i) = p;
            U = [H1 \ (D1'*p + RHS1)
                H2 \ (D2'*p + RHS2)];
            opt.U(:,i) = U;
            sol = [U;p];
        end


    end

    if study.t(i) >= next_frame_time || i==length(study.t)
        count = count + 1;
        opt.Pr(:,count) = sol(end-opt.neqnP+1:end);
        opt.Pr(:,count) = opt.Pr(:,count)-mean(opt.Pr(:,count));
        opt.U(:,count) = sol(1:2*opt.neqnV);
        opt.movie_step(count) = i;


        save_step(folderName,fname,mesh,study.N_interp,sol(1:2*opt.neqnV)'...
            ,opt.Pr(:,count),opt.neqnV,study.t(i),i,count,opt.nel)
        next_frame_time = next_frame_time + fdt;
    end
    if any(isnan(sol))

        fprintf(['STUDY FAILED -- Velocity field or Pressure field contains NaN ...' ...
            'values\nField diverged after %d iterrations out of %d'],i,ndt);
       
        return
    end


end
end

function save_step(folderName,fname,mesh,N_interp,U,Pr,neqnV,time,step,count,nel)
filename = sprintf('%s%1d.f%05d',fname,0,count);
full_path = fullfile(folderName,filename);
[data_interp] = twoD_element_interpolator(mesh,N_interp,U(1:neqnV)',U(neqnV+1:end)',Pr(:));
flag = writenek(full_path,data_interp,[N_interp,N_interp,1],1:nel,time,step,'XUP','le',4,6.54321);

if flag~=0
    fprintf('Writenek failed for timestep %d',step)
end
end

function [xDistance, yDistance] = findNearestXYDistances(points)
numPoints = size(points, 1);
xDistance = zeros(numPoints, 1);
yDistance = zeros(numPoints, 1);

% Extract unique x and y coordinates to determine grid spacing
uniqueX = unique(points(:,1));
uniqueY = unique(points(:,2));

% Calculate the minimum distances assuming uniform grid spacing
if length(uniqueX) > 1
    deltaX = min(diff(uniqueX)); % smallest x-axis spacing between neighbors
else
    deltaX = 0; % All points have the same x coordinate
end

if length(uniqueY) > 1
    deltaY = min(diff(uniqueY)); % smallest y-axis spacing between neighbors
else
    deltaY = 0; % All points have the same y coordinate
end

% Determine closest distances
for i = 1:numPoints
    currentPoint = points(i, :);

    % Find the closest point in the x direction
    if deltaX ~= 0
        xOptions = points(points(:,1) ~= currentPoint(1), :); % Exclude the current point
        [~, idxX] = min(abs(xOptions(:,1) - currentPoint(1)));
        xDistance(i) = abs(xOptions(idxX, 1) - currentPoint(1));
    end

    % Find the closest point in the y direction
    if deltaY ~= 0
        yOptions = points(points(:,2) ~= currentPoint(2), :); % Exclude the current point
        [~, idxY] = min(abs(yOptions(:,2) - currentPoint(2)));
        yDistance(i) = abs(yOptions(idxY, 2) - currentPoint(2));
    end
end
end


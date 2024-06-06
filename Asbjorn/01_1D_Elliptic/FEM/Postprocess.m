function opt = Postprocess(mesh,study,opt)
% Postprocess for stresses and strains, etc.
% function opt = Postprocess(mesh,study,opt); 
%

resolution = 1e3;
opt.elementSolution = zeros(opt.nel,resolution);
opt.elementX = zeros(opt.nel,resolution);

% Loop over elements
for e=1:opt.nel

    % Get coordinates
    x = mesh.X(mesh.IX(e,2):mesh.IX(e,end-1),2);
    
    % Get element dofs
    edof=mesh.IX(e,2):mesh.IX(e,study.N+2);

    knots = mesh.X(edof,2)';

    ydata = opt.U(edof)';

    t = linspace(x(1),x(end),resolution);

    opt.elementX(e,:) = t;
    opt.elementSolution(e,:) = LagrangeFormInterpolation(knots,ydata,t);

end

% if strcmp(study.element,'coupled')==1
%     nna = opt.neqn_a;
%     nn = opt.neqn_p;
% else
%     nn = opt.neqn/3;
% end
% 
% if strcmp(study.analysis,'static')==1
%     % Assign output parameter
%     opt.pdof = 3*mesh.wdof(1) - (3-mesh.wdof(2));
%     %p.pdof = find(p.P~=0);
%     opt.w = opt.U(opt.pdof,:);
% 
% elseif strcmp(study.analysis,'freq_direct')==1 ...
%         || strcmp(study.analysis,'freq_modal')==1 ...
%         || strcmp(study.analysis,'freq_modal_acc')==1 ...
%         || strcmp(study.analysis,'freq_modal_comp')==1 ...
%         || strcmp(study.analysis,'freq_modal_comp_acc')==1 ...
%         || strcmp(study.analysis,'freq_modal_coupled_acc')==1 ...
%         || strcmp(study.analysis,'freq_modal_coupled')==1
%     % Assign output parameter
%     opt.pdof = 3*mesh.PointLoads(1) - (3-mesh.PointLoads(2));
%     opt.centerdof = 3*floor((nn+1)/2) - (3-1);
%     %p.pdof = find(p.P~=0);
%     opt.w = opt.U(opt.pdof,:);
% 
%     opt.wtot = zeros(1,size(opt.U,2));
%     for j = 1:size(opt.U,2)
%         opt.wtot(j) = opt.U(:,j).'*opt.U(:,j);
%     end
%     opt.w_center = opt.U(opt.centerdof,:);
% 
%     try
%     specnode = 894.0; spec_dof = 3*specnode - (3-1); %Self specified point to extract data
%     opt.w_spec = opt.U(spec_dof,:);
%     catch
% 
%     end
% 
%     %----------------------------Energies:--------------------------------%
%     % Only calculate energies of study.omega has one entry
%     if length(study.omega) == 1 && ~strcmp(study.element,'acoustic')==1
% 
%         if strcmp(study.element,'mindlin')==1
%             opt.Up = opt.U;
%         elseif strcmp(study.element,'coupled')==1
%             opt.Up = opt.U(1:opt.neqn_p);
%             opt.Mp = opt.M(1:opt.neqn_p,1:opt.neqn_p);
%             opt.Kp = opt.K(1:opt.neqn_p,1:opt.neqn_p);
%             opt.Cp = opt.C(1:opt.neqn_p,1:opt.neqn_p);
%             opt.Pp = opt.P(1:opt.neqn_p);
%         end
%             opt.Up = opt.U;
%             opt.Mp = opt.M;
%             opt.Kp = opt.K;
%             opt.Cp = opt.C;
%             opt.Pp = opt.P;
%     %GLOBAL Average kinetic and potential energy
%     opt.aEkin = study.omega.^2/4.*opt.Up.'*opt.Mp*conj(opt.Up);
%     opt.aEpot = 1/4*opt.Up.'*opt.Kp*conj(opt.Up);
% 
%     %Global input and dissipated energy
%     opt.dEinp = -pi*imag(opt.Up).'*opt.Pp(:,1);
%     opt.dEdis = pi*study.omega.*opt.Up.'*opt.Cp*conj(opt.Up);
% 
%     %Local energies
%     opt.S = zeros(opt.nel,2); %CENTER - x and y
%     S_hat = zeros(opt.nel,2); %Integated - x and y
%     opt.Etra = zeros(opt.nel,2);
% 
%     for e=1:opt.nel
%         % element nodes
%         nen = mesh.IX(e,2:5);
%         % Get coordinates
%         xy = mesh.X(nen,2:3);
%         % Get element dofs
%         for i=1:4
%             edof(3*i-2) = 3*nen(i)-2;
%             edof(3*i-1) = 3*nen(i)-1;
%             edof(3*i-0) = 3*nen(i)-0;
%         end
%         %Check pressure loads
%         if isempty(mesh.PressureLoads)
%             q = 0;
%         else
%             q = mesh.PressureLoads(1,3);
%         end
% 
%         % Get material parameters
%         matID = mesh.IX(e,6);
%         if isfield(mesh,'thk')
%             thk = mesh.thk(e);
%         else
%             thk = mesh.Material(matID,1);
%         end
%         E   = mesh.Material(matID,2);
%         nu  = mesh.Material(matID,3);        rho = mesh.Material(matID,4);        
%         alpha = mesh.Material(matID,5);      beta = mesh.Material(matID,6);
%         eta = mesh.Material(matID,7);
% 
%         % Get the integrated element matrices
%         [ke,me,~] = elementMatrixMindlin(xy(:,1)',xy(:,2)',E,thk,nu,rho,study.intMethod,q);
%         ce = alpha*me + beta*ke;
% 
%         % ELEMENT ENERGIES
%         opt.aEkin_elem(e,:) = study.omega.^2/4.*opt.Up(edof).'*me*conj(opt.Up(edof));
%         opt.aEpot_elem(e,:) = 1/4*opt.Up(edof).'*ke*conj(opt.Up(edof));
%         opt.dEdis_elem(e,:) = pi*study.omega.*opt.Up(edof).'*ce*conj(opt.Up(edof));
% 
%         %POYNTING VECTOR
%         [qx,qy,Qx,Qy] = elementMatrixMindlinPoynting(xy(:,1)',xy(:,2)',E,thk,nu,study.intMethod);
%         i = sqrt(-1);
% 
%         %Evalueted at element center
%         opt.S(e,1) = 1/2*study.omega.*real(-i*(opt.Up(edof).'*qx*conj(opt.Up(edof))));
%         opt.S(e,2) = 1/2*study.omega.*real(-i*(opt.Up(edof).'*qy*conj(opt.Up(edof))));
%          if strcmp(study.element,'coupled')==1
%              opt.Sp = opt.S;
%          end
% 
%         %Intergration of the elemet
%         S_hat(e,1) = 1/2*study.omega.*real(-i*opt.Up(edof).'*Qx*conj(opt.Up(edof)));
%         S_hat(e,2) = 1/2*study.omega.*real(-i*(opt.Up(edof).'*Qy*conj(opt.Up(edof))));
% 
%         T = 2*pi/study.omega;
%         %Energy trough an element
%         opt.Etra(e,1) = S_hat(e,1)*T;
%         opt.Etra(e,2) = S_hat(e,2)*T;
%     end
%     end
% 
%         %----------------------------Energies:--------------------------------%
%     % Only calculate energies of study.omega has one entry
%     if length(study.omega) == 1 && strcmp(study.element,'coupled')==1
% 
% 
%         if strcmp(study.assembly,'standard')==1
%             % Loop over element and integrate
%             for e=1:opt.nel_a
%                 % element nodes
%                 nen = mesh.IX3D(e,2:9);
% 
%                 % Get coordinates
%                 xyz = mesh.X3D(nen,2:4);
% 
%                 % Get element dofs
%                 % simple due to only having 1-dof per node
%                 edof = nen;
% 
%                 % Get material parameters
%                 matID = mesh.IX3D(e,10);
%                 rho = mesh.Material(matID,1);        cspeed = mesh.Material(matID,2);
% 
%                 if length(mesh.Material(matID,:))>5
%                     alp = mesh.Material(matID,5);   bet = mesh.Material(matID,6);
%                 else
%                     alp = 0;   bet = 0;
%                 end
% 
%                 % Get the integrated element matrices
%                 [ke, me, fe,qx,qy,qz] = elementMatrixAcousticIntensity(xyz(:,1)',xyz(:,2)',xyz(:,3)',cspeed,rho);
% 
%                 Pres = opt.U(opt.neqn_p+1:end);
%                 Ix = -1/(2*rho*study.omega)*real(sqrt(-1)*Pres(edof).'*qx*conj(Pres(edof)));
%                 Iy = -1/(2*rho*study.omega)*real(sqrt(-1)*Pres(edof).'*qy*conj(Pres(edof)));
%                 Iz = -1/(2*rho*study.omega)*real(sqrt(-1)*Pres(edof).'*qz*conj(Pres(edof)));
% 
%                 %Poyinting vector
%                 opt.I(e,:) = [Ix Iy Iz];
% 
%             end
%         end
%     end
% 
% end
end
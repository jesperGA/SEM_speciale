clear
close all
clc
mat = [1.1,1.2,1.3,1.4;
    1.1,1.2,1.3,1.4];

GLL = 3:2:30;
% for i = 1:numel(GLL)
% n_GLL = GLL(i); %Specify number of GLL points
n_GLL = 5;

L = pi;

[xi,w,~] = lglnodes(n_GLL-1);
study.xi = xi;study.w = w;study.n_GLL = n_GLL;

[iglob, xN,yN] = MeshBox(1,1,2,2,n_GLL);
%Create mesh like  Rønquist in Figure X. 4 elements in a 2x2 constallation.
%and a sinus shaped top. 
for e = 3:4
    glob_mod = iglob(:,:,e);
    for k = 1:size(glob_mod,1)
        if e==4 && k==1
        else
            x_spec = xN(glob_mod(k,end));
            yN(glob_mod(k,end)) = 1+1/4*sin(pi*x_spec);
            ratio = (1/2+1/4*sin(pi*x_spec))/(1/2);
            yN(glob_mod(k,2:end-1)) = (yN(glob_mod(k,2:end-1))-1/2)*(ratio)+1/2;
        end
    end

end

figure()
% plotmesh(iglob,xN,yN)
mesh.IX = iglob;
mesh.X = [ones(length(xN),1),xN,yN]; %Save grid to mesh-struct
opt = [];
[opt,study] = AssemblyQuad(mesh,opt,study);

% end

% opt = [];
%
% [opt,study] = Assembly1bar(mesh,opt,study);
% xv = mesh.X(:,2);
% f = exp(xv).*(cos(xv)-sin(xv));
% A = opt.K;B = opt.M;
% f(1) = 0;f(end) = 0;
% A(1,:) = 0;A(:,1) = 0;
% A(end,:) = 0;A(:,end) = 0;
% A(1) = 1; A(end) = 1;
%
%
% u = A\(B*f);
% % figure(2)
% analytical = -exp(xv).*cos(xv)/2 - exp(xv).*sin(xv)/2 - ((1 + exp(pi)).*xv)./(2*pi) + 1/2;
% % % analytical = -sin(xv);
% err(i) = norm(analytical - u,inf);
% % plot(xv,u,xv,analytical,'LineWidth',2);
% % legend("Numerical","analytical")
%
% end
%
% figure(3);semilogy(2*GLL+1,err,'-ok','LineWidth',3)
% grid on
%
% % Setting the font size to 18
% set(gca, 'FontSize', 18);
%
% % Adding labels with LaTeX interpreter
% xlabel('$n_{dof}$', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('$\| u-u_{corr} \|_{\infty}$ Error', 'Interpreter', 'latex', 'FontSize', 18);

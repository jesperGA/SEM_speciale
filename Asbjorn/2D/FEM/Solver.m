function opt = Solver(mesh,study,opt)
% function opt = Solver(mesh,study,p)
%
% Input
% study.analysis =  'static'
%                   ! More will come from next week !
%
% study.X with X = neig (int), omega (float(:))
%
% mesh: mesh data (X,IX,bound,....)
% opt: Matrices (K,C,M,P,Null)
%
% Output
% opt.U , opt.w , (opt.pdof)


tic;
opt.U = opt.A \ (opt.f);
time = toc;
fprintf('Solution found in %f s\n',time);

end
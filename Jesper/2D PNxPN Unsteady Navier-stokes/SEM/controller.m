function [opt, study] = controller(mesh,study)

opt = [];
tic
if strcmp(study.P_order,'PnPn-2') == 1
[opt,study] = AssemblyQuad(mesh,opt,study);
elseif strcmp(study.P_order,'PnPn') == 1
    [opt,study] = PnPn_AssemblyQuad(mesh,opt,study);
end
dt = toc;
fprintf('Time for assembly: %2.3f seconds\n',dt)
if strcmp(study.study_type,'steady') == 1
    opt = solver(opt,study,mesh);
elseif strcmp(study.study_type,'unsteady') == 1
    tic
    opt = transient_solver(opt,study,mesh);
    dt = toc;
    fprintf('Time for solution: %2.3f seconds',dt)
end

end
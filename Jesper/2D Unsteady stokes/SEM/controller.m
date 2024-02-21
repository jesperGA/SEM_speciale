function [opt, study] = controller(mesh,study)

opt = [];
[opt,study] = AssemblyQuad(mesh,opt,study);
if strcmp(study.study_type,'steady') == 1
    opt = solver(opt,study,mesh);
elseif strcmp(study.study_type,'unsteady') == 1
    opt = transient_solver(opt,study,mesh);
end

end
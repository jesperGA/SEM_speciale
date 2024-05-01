function [opt, study] = controller(mesh,study)

opt = [];
para_run = 1;
try
    gcp;
catch
    para_run = 0;
    disp('Assembling in serial')
    tic
    [opt,study] = AssemblyQuad(mesh,opt,study);
    dt = toc;
end
if para_run
    tic
    [opt,study] = parAssemblyQuad(mesh,opt,study);
    dt = toc;
end
fprintf('Time for assembly: %2.3f seconds\n',dt)
if strcmp(study.study_type,'steady') == 1
    opt = solver(opt,study,mesh);
elseif strcmp(study.study_type,'unsteady') == 1
    tic
    opt = transient_solver(opt,study,mesh);
    dt = toc;
    fprintf('\nTime for solution: %2.3f seconds',dt)
end

end
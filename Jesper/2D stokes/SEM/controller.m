function [opt, study] = controller(mesh,study)

    opt = [];
    [opt,study] = AssemblyQuad(mesh,opt,study);
    opt = solver(opt,study);
    disp('Assembly done')

end
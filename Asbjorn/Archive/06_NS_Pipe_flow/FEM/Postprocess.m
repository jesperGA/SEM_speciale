function opt = Postprocess(mesh,study,opt)
% Postprocess for stresses and strains, etc.
% function opt = Postprocess(mesh,study,opt); 
fname = 'pipe';
RefineTimes = 10;
steps = 1:round(length(study.t)/100):length(study.t);
save2nek5000(fname, steps, opt,mesh,study,RefineTimes);

end
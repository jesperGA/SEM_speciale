function plotPressure(mesh,study,opt)

figure()
for p=1:opt.neqn_p
    x(p,1)=mesh.Xp(p,2);
    y(p,1)=mesh.Xp(p,3);
    z(p,1)= full(opt.p(p));
    x(p,2)=mesh.Xp(p,2);
    y(p,2)=mesh.Xp(p,3);
    z(p,2)= 0;
end
plot3(x',y',z','k')
enhance_plot(0,0,0,0,0)
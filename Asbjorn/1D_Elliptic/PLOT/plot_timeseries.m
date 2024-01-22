function plot_timeseries(mesh,opt,n_steps)
figure()
pl=plot(mesh.X(:,2),opt.u(:,1));
ylim([min(min(opt.u)) max(max(opt.u))])
xlim([0-mesh.X(end,2)*0.01 mesh.X(end,2)*1.01])
xlabel('x')
grid on
drawnow
enhance_plot(0,0,0,0,0);
for i = 1 : n_steps/1000: n_steps
    pl.YData = opt.u(:,i);
    pause(0.01)
    drawnow
end
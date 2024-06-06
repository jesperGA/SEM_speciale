R = 69; umax = 1.5;
r = linspace(0,2*R,1e5);

u = umax.*(1-((r-R)/R).^2);
mean(u)
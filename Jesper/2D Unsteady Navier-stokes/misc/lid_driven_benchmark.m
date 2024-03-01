function [errorV] = lid_driven_benchmark(y_in,u_in)

A = readmatrix("misc\lid_benchmark.txt");

%Load benchmark results
benchmark2 = A;
y_db =  benchmark2(:,1);
% px_bench = benchmark2(:,6) ;
p_db = benchmark2(:,3);
% v_bench = benchmark2(:,5);
u_db = benchmark2(:,2);

% p = interp1(y_in,p_in,y_db);
u = interp1(y_in,u_in,y_db);

% errorP = norm(p-p_db,inf);
errorV = norm(u-u_db,inf);
end

fprintf('Solver status:') %For AeStheTics
prev_mes = 0;
ndt = 100000;
for i = 2:ndt
    %% For AeStheTics
    if mod(i,500) ==0
        status_procent = (i/ndt)*100;
        mess = sprintf('%3.0f',status_procent);
        back_space = repmat('\b',1,prev_mes);
        fprintf(back_space)
        fprintf(mess)
        prev_mes = length(mess);
        pause(0.1)
    end
    %%

end
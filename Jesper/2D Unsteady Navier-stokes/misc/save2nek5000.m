function save2nek5000(fname, T, opt,mesh,study,N_interp)
% save2nek5000(fname, T, opt,mesh,N_interp)
% Saves visualisation to the timesteps in the vector T (ie T = 1:10 for the first 10 timesteps)
% The solution is intepreted to N_interp X N_interp points in each element.
% Visualisation saved in the folder fname_anis\*
fprintf('\nSAVING TO ANIMATIONS of %s',fname)
N0 = 1;
Nt = numel(T);

% Create a directory for the files
folderName = [fname, '_anis'];  % Construct folder name
if ~exist(folderName, 'dir')
    mkdir(folderName);          % Create the folder if it doesn't exist
end

file_path = fullfile(folderName, [fname, '.nek5000']);

% Open the file for writing with .nek5000 extension
fileID = fopen(file_path, 'w+');

% Write the specified lines with the variable contents into the file
fprintf(fileID, 'filetemplate: %s%%01d.f%%05d\n', fname);
fprintf(fileID, 'firsttimestep: %d\n', N0);
fprintf(fileID, 'numtimesteps: %d\n', Nt-1);

% Close the file
fclose(fileID);
counter = 0;
counter(T) = 1:numel(T);

t = study.t; U = opt.U;Pr = opt.Pr; neqnV = opt.neqnV;nel = opt.nel;
for k = 1:numel(T)
    i = T(k);
    count = counter(i);
    time = t(i);
    filename = sprintf('%s%1d.f%05d',fname,0,count);
    full_path = fullfile(folderName,filename);
    [data_interp] = twoD_element_interpolator(mesh,N_interp,U(1:neqnV,k),U(neqnV+1:end,k),Pr(:,k));
    flag = writenek(full_path,data_interp,[N_interp,N_interp,1],1:nel,time,i,'XUP','le',4,6.54321);

    if flag~=0
        fprintf('Writenek failed for timestep %d',i)
    end
end
fprintf('\n Done saving animation')
end
function save2nek5000(fname, steps, opt,mesh,study,RefineTimes)
% save2nek5000(fname, T, opt,mesh,N_interp)
% Saves visualisation to the timesteps in the vector T (ie T = 1:10 for the first 10 timesteps)
% The solution is intepreted to N_interp X N_interp points in each element.
% Visualisation saved in the folder fname_anis\*

N0 = 1;
Nt = numel(steps);
N_interp = (study.N+1)*RefineTimes;
nel = opt.nel;

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
t = study.t;
IX = mesh.IX;
IXp = mesh.IXp;
X = mesh.X;
Xp = mesh.Xp;
N = study.N;
u1 = opt.u1;
u2 = opt.u2;
p = opt.p;

parfor i = 1:length(steps)
    time = t(steps(i));
    filename = sprintf('%s_anis\\%s%1d.f%05d',fname,fname,0,i);
    [data_interp] = twoD_element_interpolator(N_interp,nel,u1(:,steps(i)),u2(:,steps(i)),p(:,steps(i)), IX, IXp, X, Xp, N);
    flag = writenek(filename,data_interp,[N_interp,N_interp,1],1:nel,time,i,'XUP','le',4,6.54321);

    if flag~=0
        fprintf('Writenek failed for timestep %d',i)
    end
end
end
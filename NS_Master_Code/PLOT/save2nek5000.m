function save2nek5000(fname, opt, mesh, study, RefineTimes)
% save2nek5000(fname, opt, mesh, study, RefineTimes)
% Saves visualization to the timesteps in the vector t specified by t_steps_coarse
% Visualization saved in the folder specified by study.folderName

t = study.t;
IX = mesh.IX;
IXp = mesh.IXp;
X = mesh.X;
Xp = mesh.Xp;
N = study.N;
u1 = opt.u1;
u2 = opt.u2;
p = opt.p;
t_steps = study.t_steps;
t_steps_coarse = 1:round(length(t_steps)/100):length(t_steps);

N0 = 1;
Nt = numel(t_steps_coarse);
N_interp = (study.N+1) * RefineTimes;
nel = opt.nel;

% Ensure the base folder exists
if ~exist(study.folderName, 'dir')
    mkdir(study.folderName);
end

% Create a subdirectory for the files
subFolderName = fullfile(study.folderName, fname);
if ~exist(subFolderName, 'dir')
    mkdir(subFolderName);
end

file_path = fullfile(subFolderName, [fname, '.nek5000']);

% Open the file for writing with .nek5000 extension
fileID = fopen(file_path, 'w+');

% Write the specified lines with the variable contents into the file
fprintf(fileID, 'filetemplate: %s%%01d.f%%05d\n', fullfile(fname));
fprintf(fileID, 'firsttimestep: %d\n', N0);
fprintf(fileID, 'numtimesteps: %d\n', Nt - 1);

% Close the file
fclose(fileID);

% Write data files in parallel
parfor i = 1:length(t_steps_coarse)
    time = t(t_steps(t_steps_coarse(i)));
    filename = fullfile(subFolderName, sprintf('%s%1d.f%05d', fname, 0, i));
    [data_interp] = twoD_element_interpolator(N_interp, nel, u1(:, t_steps_coarse(i)), u2(:, t_steps_coarse(i)), p(:, t_steps_coarse(i)), IX, IXp, X, Xp, N);
    flag = writenek(filename, data_interp, [N_interp, N_interp, 1], 1:nel, time, i, 'XUP', 'le', 4, 6.54321);

    if flag ~= 0
        fprintf('Writenek failed for timestep %d', i);
    end
end
end

% Define the directory where your PNG files are stored
folderPath = 'C:\Users\jga91\OneDrive - Danmarks Tekniske Universitet\SEM Speciale\Kode\SEM_speciale\Jesper\2D Unsteady Navier-stokes\misc\pipe_whole_20x4_gll5';  % Update this to your folder path

% Get a list of all PNG files in the directory
files = dir(fullfile(folderPath, '*.png'));
numFiles = length(files);

% Extract numbers from filenames for sorting
fileNumbers = zeros(numFiles, 1);
for i = 1:numFiles
    % Assuming filenames contain digits that can be sorted numerically
    digits = regexp(files(i).name, '\d+', 'match');
    if ~isempty(digits)
        fileNumbers(i) = str2double(digits{1});
    end
end

% Sort files based on extracted numbers
[~, sortOrder] = sort(fileNumbers);
sortedFiles = {files(sortOrder).name};

% Define the name of the output GIF file
outputGIF = 'pipe_Whole_20x4_GLL5.gif';

% Loop through each file, read the image, and write to GIF
for k = 1:numFiles
    % Full path to the image file
    filename = fullfile(folderPath, sortedFiles{k});
    
    % Read image
    img = imread(filename);
    
    % Capture the frame for the GIF
    [A, map] = rgb2ind(img, 256);
    
    % Write frame to GIF file
    if k == 1
        imwrite(A, map, outputGIF, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, outputGIF, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

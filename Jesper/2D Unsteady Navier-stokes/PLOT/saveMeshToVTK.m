function saveMeshToVTK(x, y, filename)
    % This function saves mesh data to a VTK file for visualization in ParaView
    % Inputs:
    %   x, y   - Mesh coordinates
    %   filename - Name of the VTK file to create

    % Open file
    fileID = fopen(filename, 'w');
    
    % Write VTK header and grid data type
    fprintf(fileID, '# vtk DataFile Version 3.0\n');
    fprintf(fileID, '2D Mesh Data\n');
    fprintf(fileID, 'ASCII\n');
    fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
    fprintf(fileID, 'DIMENSIONS %d %d 1\n', length(x), length(y));
    fprintf(fileID, 'POINTS %d float\n', length(x) * length(y));
    
    % Assuming a simple grid (Nx by Ny)
    Nx = length(x);
    Ny = length(y);
    for j = 1:Ny
        for i = 1:Nx
            fprintf(fileID, '%f %f 0.0\n', x(i), y(j));
        end
    end
    
    % Close the file
    fclose(fileID);
end
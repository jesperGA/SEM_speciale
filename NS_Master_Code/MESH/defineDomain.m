function [LX, LY, NELX, NELY] = defineDomain(study)
    % defineDomain Define the physical domain and mesh parameters based on the study example
    % Inputs:
    %   study - Struct containing the example type
    % Outputs:
    %   LX - Length of the domain in the x-direction
    %   LY - Length of the domain in the y-direction
    %   NELX - Number of elements in the x-direction
    %   NELY - Number of elements in the y-direction

    example = study.example;

    % Parameters of the square box domain and mesh
    switch example
        case {'TaylorVortex', 'LidDriven'}
            LX = 2; % Side-length of the box
            NELX = 2; % Number of elements along each side
            LY = LX;
            NELY = NELX;
        case 'Duct'
            includeXOutOf10Elem = 10;
            LX = 2.05 * includeXOutOf10Elem / 10;
            LY = 0.41;
            NELY = 10;
            NELX = 5 * NELY;
        otherwise
            error('Unknown example type');
    end

    % Ensure square elements
    if LX / NELX ~= LY / NELY
        error('Non-square elements not properly implemented. Check the C or D operator...');
    end
end

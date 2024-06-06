function [LX, LY, NELX, NELY] = defineDomain(study,polOrder,domain)

    example = study.example;
    % Parameters of the square box domain and mesh
    if strcmp(example,'Roenquist_NS') || strcmp(example,'LidDriven')
        LX = 2; % Side-length of the box
        NELX = 2; % Number of elements along each side
        LY = LX;
        NELY = LX;
    elseif  strcmp(example,'Pipe')

        includeXOutOf10Elem = 10;

        LX = 2.05*includeXOutOf10Elem/10;
        LY =0.41;
        
        nNodesX = domain * 10 * polOrder + 1;
        nNodesY = domain * 2  * polOrder + 1;
        if strcmp(study.FEM,'LFEM')
            NELX = (nNodesX - 1) * includeXOutOf10Elem/10; % Number of elements along each side
            NELY = (nNodesY - 1);
        elseif strcmp(study.FEM,'QFEM')
            NELX = (nNodesX - 1)/2 * includeXOutOf10Elem/10; % Number of elements along each side
            NELY = (nNodesY - 1)/2;
        elseif strcmp(study.FEM,'SEM')
            NELX = domain * 10 * includeXOutOf10Elem/10; % Number of elements along each side
            NELY = domain * 2;
        end

    end

end
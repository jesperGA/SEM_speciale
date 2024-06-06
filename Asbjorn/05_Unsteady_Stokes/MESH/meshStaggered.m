function [mesh] = meshStaggered(study, LX, LY, NELX, NELY, N);

    distribution=0;
    [mesh] = mesh2D(LX, LY, NELX, NELY, N, distribution);

    distribution=2;
    [meshP] = mesh2D(LX, LY, NELX, NELY, N-2, distribution);

    if strcmp(study.example,'Roenquist')==1 || strcmp(study.example,'Roenquist_Poisson')==1
        mesh.X(:,2:3)=mesh.X(:,2:3)-LX/2;
        meshP.X(:,2:3)=meshP.X(:,2:3)-LX/2;
    end

    mesh.Xp = meshP.X;
    mesh.IXp = meshP.IX;

    %-------------------------------------------------------------------------%
    %                             Material parameters                         %
    %-------------------------------------------------------------------------%
    
    if study.unsteady
        mesh.Material = [1, 1.2, 1.3, 1.4];
    else
        mesh.Material = [1, 1.2, 0, 1.4];
    end

end
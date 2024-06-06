function [mesh] = meshStaggered(study, LX, LY, NELX, NELY, N);

    distribution=0;
    [mesh] = mesh2D(LX, LY, NELX, NELY, N, distribution);

    distribution=0;
    [meshP] = mesh2D(LX, LY, NELX, NELY, N, distribution);

    if strcmp(study.example,'Roenquist') || strcmp(study.example,'Roenquist_Poisson') || strcmp(study.example,'Roenquist_NS')
        mesh.X(:,2:3)=mesh.X(:,2:3)-LX/2;
        meshP.X(:,2:3)=meshP.X(:,2:3)-LX/2;
    end

    mesh.Xp = meshP.X;
    mesh.IXp = meshP.IX;

    mesh.Material = [1, 1];

end
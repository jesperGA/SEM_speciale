function plotSolContour(mesh,study,opt,resolution,NELX, NELY)
    fntsize=22;
    N = study.N;

    name='$u_1$';
    elementLoop(opt.u1, mesh.X, mesh.IX, N, opt.nel, resolution, name, fntsize,NELX, NELY,1)
    saveas(gcf,['contourSol_',name],'epsc')
    name='$u_2$';
    elementLoop(opt.u2, mesh.X, mesh.IX, N, opt.nel, resolution, name, fntsize,NELX, NELY,1)
    saveas(gcf,['contourSol_',name],'epsc')
    name='$p$';
    elementLoop(opt.p, mesh.Xp, mesh.IXp, N-2, opt.nel, resolution, name, fntsize,NELX, NELY,2)
    xlim([min(mesh.X(:,2)) max(mesh.X(:,2))])
    ylim([min(mesh.X(:,3)) max(mesh.X(:,3))])
    saveas(gcf,['contourSol_',name],'epsc')
end

function elementLoop(sol, meshX, meshIX, N, nel,resolution, name, fntsize,NELX, NELY, distribution)
    fig = figure();
    hold on

    if distribution==1
        xx = zeros(NELX*N+1, NELY*N+1);
        yy = zeros(NELX*N+1, NELY*N+1);
        uu = zeros(NELX*N+1, NELY*N+1);
    elseif distribution==2
        xx = zeros(NELX*(N+1), NELY*(N+1));
        yy = zeros(NELX*(N+1), NELY*(N+1));
        uu = zeros(NELX*(N+1), NELY*(N+1));
    end

    e=0;
    for ey=1:NELY
        for ex=1:NELX
            
            e = e+1;
            % Get element data
            [xy, edof] = getElementData(meshX, meshIX, e);
    
            x = reshape(xy(:, 1), N + 1, N + 1);
            y = reshape(xy(:, 2), N + 1, N + 1);
            u = reshape(sol(edof), N + 1, N + 1);
    
            if distribution==1
                xx((ex-1)*N+1:ex*N+1,(ey-1)*N+1:ey*N+1) = x;
                yy((ex-1)*N+1:ex*N+1,(ey-1)*N+1:ey*N+1) = y;
                uu((ex-1)*N+1:ex*N+1,(ey-1)*N+1:ey*N+1) = u;
            elseif distribution==2
                xx((ex-1)*(N+1)+1:ex*(N+1),(ey-1)*(N+1)+1:ey*(N+1)) = x;
                yy((ex-1)*(N+1)+1:ex*(N+1),(ey-1)*(N+1)+1:ey*(N+1)) = y;
                uu((ex-1)*(N+1)+1:ex*(N+1),(ey-1)*(N+1)+1:ey*(N+1)) = u;
            end
    
        end
    end


    contourf(xx,yy,uu,resolution, 'LineColor', 'none')
    colorbar, axis image;
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fntsize);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fntsize);
    enhance_plot(0,fntsize,0,0,0)
end

function [xy, edof] = getElementData(X, IX, e)
    nen = IX(:,:,e);
    xy = X(nen, 2:3);
    edof = reshape(nen, [], 1);
end
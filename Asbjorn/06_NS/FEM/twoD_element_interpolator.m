function [data] = twoD_element_interpolator(N_interp,nel,u1,u2,p, IX, IXp, X, Xp, N);

data = zeros(nel, (N_interp)^2, 5);

XeI  = zeros(N_interp, N_interp, nel);
YeI  = zeros(N_interp, N_interp, nel);
u1eI = zeros(N_interp, N_interp, nel);
u2eI = zeros(N_interp, N_interp, nel);
peI =  zeros(N_interp, N_interp, nel);

out1 = zeros(nel,(N_interp)^2);
out2 = zeros(nel,(N_interp)^2);
out3 = zeros(nel,(N_interp)^2);
out4 = zeros(nel,(N_interp)^2);
out5 = zeros(nel,(N_interp)^2);

for e = 1:nel
    xy = X(IX(:,:,e), 2:3);
    xyp = Xp(IXp(:,:,e), 2:3);

    XX = reshape(xy(:, 1), N + 1, N + 1);
    YY = reshape(xy(:, 2), N + 1, N + 1);

    if N==1 || N==2
        XXp = reshape(xyp(:, 1), N - 1 + 1, N - 1 + 1);
        YYp = reshape(xyp(:, 2), N - 1 + 1, N - 1 + 1);    
    else
        XXp = reshape(xyp(:, 1), N - 2 + 1, N - 2 + 1);
        YYp = reshape(xyp(:, 2), N - 2 + 1, N - 2 + 1);
    end

    xxI = linspace(min(xy(:, 1)),max(xy(:, 1)),N_interp);
    yyI = linspace(min(xy(:, 2)),max(xy(:, 2)),N_interp);
    [XeI(:,:,e),YeI(:,:,e)] = meshgrid(xxI,yyI);

    [u1eI(:,:,e)] = lagrange2D(XX', YY', N, reshape(u1(IX(:,:,e)),N+1,N+1), XeI(:,:,e), YeI(:,:,e));
    [u2eI(:,:,e)] = lagrange2D(XX', YY', N, reshape(u2(IX(:,:,e)),N+1,N+1), XeI(:,:,e), YeI(:,:,e));
    if N==1 || N==2
        [peI(:,:,e)]  = lagrange2D(XXp',YYp',N-2, reshape(p(IXp(:,:,e)),N-1+1,N-1+1),XeI(:,:,e), YeI(:,:,e));    
    else
        [peI(:,:,e)]  = lagrange2D(XXp',YYp',N-2, reshape(p(IXp(:,:,e)),N-2+1,N-2+1),XeI(:,:,e), YeI(:,:,e));
    end
    
    out1(e,:) = reshape(XeI(:,:,e),1,numel(XeI(:,:,e)));
    out2(e,:) = reshape(YeI(:,:,e),1,numel(YeI(:,:,e)));
    out3(e,:) = reshape(u1eI(:,:,e),1,numel(u1eI(:,:,e)));
    out4(e,:) = reshape(u2eI(:,:,e),1,numel(u2eI(:,:,e)));
    out5(e,:) = reshape(peI(:,:,e),1,numel(peI(:,:,e)));

    data(e,:,1) = out1(e,:);
    data(e,:,2) = out2(e,:);
    data(e,:,3) = out3(e,:);
    data(e,:,4) = out4(e,:);
    data(e,:,5) = out5(e,:);
end

end
function [A,B,fe] = elementMatrixBarFEM(x);

    fe=zeros(2,1);

    % Initialize element matrices
    A = zeros(2, 2);
    B = zeros(2, 2);

    % Integrate mass matrix
    meGauss = 2;
    [megaussPoints,mew] = gaussFunc(meGauss);
    for i = 1:meGauss
        xi = megaussPoints(i);
        [Bm,J,N] = shapeFuncBar(x,xi);
        % Calculate the element stiffness matrix (ke) and mass matrix (me) numerically
        A = A + exp(x(1)+(xi+1)*J) * (Bm' * Bm) * J * mew(i);
        B = B + (N' * N) * J * mew(i);
    end
end

% Shapefunction and strain displacement for a quadrilateral
function [Bm,J,N] = shapeFuncBar(x,xi)

    L = x(2)-x(1);

    % Zero the strain-displacement, jacobian, shapefunction matrices
    Bm = zeros(1,2);
    N = zeros(1,2);
    
    N = [(1 - xi) / (2), (1 + xi) / (2)];
    Bm = [-1/L, 1/L];
    J = L / 2;

end

% First three sets of GP and W for a quadrilateral
function [gaussVec,w] = gaussFunc(gaussPoints)
    Gauss = [0 0 0 ;
            -1/sqrt(3) 1/sqrt(3) 0;
            -sqrt(0.6) 0 sqrt(0.6)];

     ww = [2 0 0;
           1 1 0;
           5/9 8/9 5/9];

    gaussVec = Gauss(gaussPoints,1:gaussPoints);
    w = ww(gaussPoints,1:gaussPoints);
end
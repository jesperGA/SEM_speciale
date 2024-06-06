function study = createStudy(study)

    if strcmp(study.example,'Roenquist')
        % Ronquist example
        study.U1 = @(X1,X2) 1-X2.^2;
        study.U2 = @(X1,X2) X2.*0;
        study.P  = @(X1,X2) sinpi(X1).*sinpi(X2);
        study.F1 = @(X1,X2) 2 + pi*cospi(X1).*sinpi(X2);
        study.F2 = @(X1,X2)     pi*sinpi(X1).*cospi(X2);
    elseif strcmp(study.example,'Roenquist_Poisson')
        % Ronquist Poisson example
        study.U1 = @(X1,X2) sin(X1) .* exp(-X2);
        study.U2 = @(X1,X2) X1 .* 0;
        study.P  = @(X1,X2) X1 .* 0;
        study.F1 = @(X1,X2) X1 .* 0;
        study.F2 = @(X1,X2) X2 .* 0;
    elseif strcmp(study.example,'Bercovier_1')
        % Bercovier (1)
        study.U1 = @(X1,X2) -256 .* X1 .^2 .* (X1 - 1) .^ 2 .* X2 .* (X2 - 1) .*(2 .* X2 - 1);
        study.U2 = @(X1,X2) study.U1(X2,X1); % Wrong sign!?
        study.P  = @(X1,X2) X1 .* 0; 
        study.F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
                            2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2));
        study.F2 = @(X1,X2) study.F1(X2,X1);
    elseif strcmp(study.example,'Bercovier_2')
        % Bercovier (2)
        study.U1 = @(X1,X2) -256 .* X1 .^2 .* (X1 - 1) .^ 2 .* X2 .* (X2 - 1) .*(2 .* X2 - 1);
        study.U2 = @(X1,X2) -study.U1(X2,X1); 
        study.P  = @(X1,X2) (X1-1/2) .* (X2-1/2);
        study.F1 = @(X1,X2) 128 .* (X1.^2 .* (X1 - 1) .^2 .*12 .* (2 .* X2 - 1) + ...
                            2 .* (X2 - 1) .* (2 .* X2 - 1) .* X2 .* (12 * X1 .^ 2 - 12 .* X1 + 2)) + X2 - 1/2;
        study.F2 = @(X1,X2) study.F1(X2,X1);
    elseif strcmp(study.example,'LidDriven')
        % Bercovier (2)
        study.U1 = @(X1,X2,t) X2 == 1;
        study.U2 = @(X1,X2,t) X1 .* 0;
        study.P  = @(X1,X2,t) X1 .* 0;
        study.F1 = @(X1,X2,t) X1 .* 0;
        study.F2 = @(X1,X2,t) X2 .* 0;
    elseif strcmp(study.example,'Roenquist_NS')
        % Ronquist example
        study.U1 = @(X1,X2,t) -cos(X1).*sin(X2).*exp(-2*t);
        study.U2 = @(X1,X2,t)  sin(X1).*cos(X2).*exp(-2*t);
        study.P  = @(X1,X2,t) -1/4 * (cos(2*X1) + cos(2*X2)).*exp(-4*t);
        study.F1 = @(X1,X2) X1 .* 0;
        study.F2 = @(X1,X2) X1 .* 0;
    end
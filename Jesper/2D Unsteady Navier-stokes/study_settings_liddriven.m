function study = study_settings_liddriven

study.solve_type = 'direct'; %uzawa
if strcmp(study.solve_type,'direct') == 1
    study.direct_type = 'backslash';
    % study.direct_type = 'LU';
end
% study.solve_type = 'uzawa';
study.study_type = 'unsteady';
% study.precon = 'mhat';
study.precon = 'P';
% study.study_t ype = 'steady';

if strcmp(study.study_type,'unsteady') == 1
    study.T = 0.5;
    % study.nt = 1e3;
    % study.t = linspace(0,study.T,study.nt);
    study.t = 0:1e-5:study.T;
    study.nt = length(study.t);
    study.dt = (study.t(2)-study.t(1));

    % study.int_type = 'BDFk'; %Equivalent of solving Unsteady stokes.
    % study.int_type = 'BDF3EX3'; %First order bdf for linear terms. 3 order for nonlinear terms.
    study.int_type = 'explicit';
    % study.RE = 100;
    % study.BDF_order = 1;

    study.U10 = 0;
    study.U20 = 0;

    % study.BC_type = 'dynamic';
    study.BC_type = 'static';
    study.p_type = 'liddriven';
    % study.p_type = 'ronquist';

    study.Brinkmann = 0;
    study.alpha = 1e6;
    % study.fps = round(1/study.dt);
    study.fps = 200;
end


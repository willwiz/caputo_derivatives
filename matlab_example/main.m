%%  Define some parameters ...
    clear

%     polyp  = [7.263437  -19.71696  100.7997  -509.0422  1106.240  -1204.629  718.3882  -236.1940  39.39095  -2.497073];
    polyp  = [0.5];
    par    = [1.0];
    period = 1;
    delta  = 10;
    alp    = 10;
    de     = 2;
    np     = 9;

%%  Load test kinematics ...
    import kinematics.*
    test = load_poly_sin_2D_x([0.5, 0.0, 0.0, -0.5], polyp, 2 );

%%  Import matlaws ...
    import matlaws.*
    import frac.*

%%   Run simulation
    for de = [1,2,3]
        % paramaters of the simulation
        Tf   = 2.0;
        nt   = 2*10^de;
        dt   = Tf/nt;
        k    = 1; %max(floor(nt/20),1);
        a    = alp/100.0;
        d    = delta/100.0;
        tag  = strcat('polyx_',string(period),'_',string(delta),'_',string(alp),'_',string(np),'_',string(de));
        % run the simulation
        t = cputime;
        disp(['Running ' test.name])
        [time, defgrad, fval, Dfval, storeR, storeL] = diffeq_approx(@law_LE_short, par, test, d, a, 9, Tf, dt, k);
        t2 = cputime;
        disp(['  - finished running in ' num2str(t2-t)])

        %   Post Processing
        disp(['Post Processing ' test.name ' ' tag])
        if ~exist(fullfile('output',tag), 'dir')
           mkdir(fullfile('output',tag))
        end
        [stress, pres] = tool_hydrostatic_pressure_2D(Dfval, 3, 3);
        %   Export Stuff
        dlmwrite(fullfile('output',tag,'vF.txt'),     defgrad,  'delimiter','\t','precision', '%18.12f')
        dlmwrite(fullfile('output',tag,'time.txt'),   time,     'delimiter','\t','precision', '%18.12f')
        dlmwrite(fullfile('output',tag,'fval.txt'),   fval,     'delimiter','\t','precision', '%18.12f')
        dlmwrite(fullfile('output',tag,'stress.txt'), stress,   'delimiter','\t','precision', '%18.12f')
        disp(['  - finished writing in ' num2str(cputime-t2)])
    end






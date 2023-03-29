function [time, x, fn, Df, storeR, storeL] = diffeq_approx(f, par, arg, delta, alp, N, T_f, dt, kn)
  % This function simulates the caputo differential equation for the function f on [0,T_f]
  % using the method defined in fracp and the time-step size dt.
  %
  %    f - a function we are to take the fractional derivative of, e.g.
  %        f = f(t, arg), where t is the time 
  %
  %    par - parameters of the model Sv, such as modulus
  %
  %    arg - e.g. structure describing the the strain/deformation gradient
  %    see +kinematics
  %
  %    delta - scaling paramter on LHS
  %
  %    alpha - fractional derivative order
  %
  %    N - number of prony terms
  %
  %    T_f - the final time of our simulation ...
  %
  %    dt - the time-step size ...
  %
  %    fracp - the structure that stores how to do fractional derivative
  %            approximation as well as the test structure ...
  %
  %    kn - how often to export results
  %
  % by David Nordsletten, 2018
  % updated by Will Zhang, 2020
  %
    % Figure out how many iterations to simulate
    N_steps  = max(round(T_f / dt),1);
    N_export = floor(N_steps/kn) + 1;
    % Initialize Caputo Parameters
    fracpL = frac.caputo_approx_initialize(alp, T_f, N);
    fracpR = frac.caputo_approx_initialize(alp, T_f, N);
    % Determine dimension of variables
    x0 = arg.arg_func(0);
    fracpR.f_prev = f(par, x0);
    fracpL.f_prev = f(par, x0);
    x0 = reshape(x0',1,[]);
    x_size = length(x0);
    f_size = length(fracpR.f_prev);
    % Set up output arrays
    time  = zeros(N_export, 1);
    x     = zeros(N_export, x_size);
    fn    = zeros(N_export, f_size);
    Df    = zeros(N_export, f_size);
    fracpR.Q = zeros(fracpR.N, f_size);
    fracpL.Q = zeros(fracpL.N, f_size);
    storeR = zeros(N_export, (N + 1)*f_size);
    storeL = zeros(N_export, (N + 1)*f_size);
    fracpL.delta = delta;
    % Set up the initial conditions
    t0      = 0.0;
    time(1) = t0;
    x(1,:)  = x0;
    fn(1,:) = fracpL.f_prev;
    Df(1,:) = fracpL.f_prev;
    storeR(1,:) = [fracpR.f_prev reshape(fracpR.Q',1,[])];
    storeL(1,:) = [fracpL.f_prev reshape(fracpL.Q',1,[])];
    % Set up export counter
    k = 1;
    for i = 1:N_steps % go through the main time loop ...
        t0 = dt + t0; % update the current time ...
        x1 = arg.arg_func(t0); % compute the arg, i.e. the strain
        % compute the right hand side and update
        [v, fracpR] = frac.caputo_approx_iter(f, par, x1, fracpR, dt);
        % solve the left hand side, i.e. S_star in the paper
        [v, fracpL] = frac.diffeq_approx1_iter(v, fracpL, dt);
        % Filling out the export data arrays ...
        if rem(i, kn)  == 0
            k = k + 1;
            time(k,1)  = t0;
            x(k,:)     = reshape(x1',1,[]);
            fn(k,:)    = fracpR.f_prev;
            Df(k,:)    = v;
            storeR(k,:) = [fracpR.f_prev reshape(fracpR.Q',1,[])];
            storeL(k,:) = [fracpL.f_prev reshape(fracpL.Q',1,[])];
        end
        
    end
end
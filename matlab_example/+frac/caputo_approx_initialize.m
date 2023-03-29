function [fracp] = caputo_approx_initialize(alp, Tf, N)
  % Imports the Prony Serie parameter for the problem given
  %       alpha
  %       TF : time interval of the simulation
  %       N_p : Number of prony terms
  % The scaling of the parameters were described in the paper
  load(['coeffs-opt-refined-100steps-' num2str(N) '-500.mat'], 'betam', 'taum');
  cycle_freq = 2.0 * pi / Tf;
  dx_alp = 1 / size(betam,2);
  alpm = dx_alp:dx_alp:1;
  tau  = zeros(N,1);
  beta = zeros(N,1);
  for i = 1:N
     tau(i,1) = interp1(alpm,taum(i,:),alp) / cycle_freq;
     beta(i,1) = interp1(alpm,betam(i,:),alp) * cycle_freq^alp;
  end
  beta(N+1,1) = interp1(alpm,betam(N+1,:),alp) * cycle_freq^(alp-1);

  % Building necessary approximation quantities ...
  fracp.N = N;
  fracp.alp = alp;
  fracp.Tf  = Tf;
  fracp.tau = tau;
  fracp.beta = beta;
  fracp.Q = [];

end
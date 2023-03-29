function [ v, fracp ] = caputo_approx_iter(f, par, arg, fracp, dt)
  % This function computes the update for 1 iteration of the caputo
  % derivative of the function f on [0,T_f] with the step size of dt
  %
  %    f - a function we are to take the fractional derivative of, e.g.
  %        f = f(t, fracp), where t is the time and fracp is a struc storing details of the function and scheme ...
  %
  %    par - parameters of the model Sv, such as modulus
  %
  %    arg - e.g. the strain/deformation gradient
  %
  %    fracp - the structure that stores how to do fractional derivative
  %            approximation as well as the test structure ...
  %
  %    dt - time step size
  %
  % by David Nordsletten, 2018
  % updated by Will Zhang, 2020
  %
  % determine the size of our function f ...
  fn   = f(par, arg);
  df   = (fn - fracp.f_prev);
  dfdt = df / dt;
  v    = fracp.beta(end) * dfdt;
  for j = 1:fracp.N % Looping over the number of prony terms ...
  % Commented section implements the midpoint rule rather than backward
  % Euler
%     e = exp(- 0.5 * dt / fracp.tau(j)); % calculate the scale factor ...
%     fracp.Q(j,:) = e * e * fracp.Q(:,j) + e * fracp.beta(j) * (fn - fnm); % updating the prony loss terms ...
    fracp.Q(j,:) = fracp.tau(j) / (fracp.tau(j) + dt) * (fracp.Q(j,:) + fracp.beta(j) * df); % updating the prony loss terms ...
    v = v + fracp.Q(j,:); % updating the fractional derivative approximation ...
  end
  fracp.f_prev = fn;
end
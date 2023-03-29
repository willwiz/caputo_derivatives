function [ v, fracp ] = diffeq_approx1_iter(fn, fracp, dt)
  % This function computes and updates the caputo differential equation for
  % the function f for 1 iteration using the method defined in fracp and 
  % the time-step size dt.
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
  %.
  % Setup
  coeff = zeros(fracp.N);
  c0 = fracp.beta(end) / dt;
  for j = 1:fracp.N % Looping over the number of prony terms ...
    coeff(j) = fracp.tau(j) / (fracp.tau(j) + dt);
    c0 = c0 + fracp.beta(j)*coeff(j);
  end
  c0 = fracp.delta*c0;
  % Do solution
  v = fn;
  for j = 1:fracp.N
      v = v - fracp.delta * coeff(j) * fracp.Q(j,:);
  end
  v = (v + c0 * fracp.f_prev)/(1.0 + c0);
  % Do some updates for stores
  df           = v - fracp.f_prev;
  fracp.f_prev = v;
  for j = 1:fracp.N
    fracp.Q(j,:) = coeff(j) * (fracp.Q(j,:) + fracp.beta(j) * df); % updating the prony loss terms ...
  end
end
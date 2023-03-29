function [stress] = law_LE(par, F)
   stress = par(1) * (F + F' - 2.0*eye(2));
   stress = reshape(stress', 1, []);
end
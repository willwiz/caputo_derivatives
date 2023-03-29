function [stress] = law_LE_short(par, F)
   stress = par(1) * (F + F' - 2.0*eye(2));
   stress = [stress(1,1) stress(1,2) stress(2,2)];
end
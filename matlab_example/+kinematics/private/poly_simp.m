function f = poly_simp(pars, t)
    tk = t;
    f  = pars(1)*tk;
    for i = 2:length(pars)
        tk = tk * t;
        f = f + pars(i)*tk;
    end
end
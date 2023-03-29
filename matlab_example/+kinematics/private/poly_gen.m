function f = poly_gen(pars, expo, t)
    if (size(pars)~=size(expo)) 
        error("parameters and exponents have different sizes") 
    end
    f = pars(1)*t^expo(1);
    for i = 2:size(pars)
        f = f + pars(i)*t^expo(i);
    end
end
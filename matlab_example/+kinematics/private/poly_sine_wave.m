function f = poly_sine_wave(a, n, t)
    t2 = -t*t;
    tk = t;
    f  = tk;
    if n > 1
        for k = 2:n
            tk = tk*t2;
            f  = f + tk/factorial(2*k-1);
        end
    end
    f = a*f;
end


# Szymon Pilawka 254649

module RootApprox

    #= = =
    # funkcja znajduje zero funkcji metodą bisekcji (podziału)
    #
    # f - funkcja
    # a - początek przedziału
    # b - koniec przedziału
    # delta - dokładność obliczeń - przedział
    # epsilon - dokładność obliczeń - wartość
    #
    # return (r: przybliżenie pierwiastka, v: wartość w r, it: ilość iteracji, err: flaga błędu)
    = = =#
    function mbisekcji(f, a::Float64, b::Float64, delta::Float64, epsilon::Float64)
        # wartość funkcji w punkcie a
        local fa = f(a)
        # wartość funkcji w punkcie b
        local fb = f(b)
        # ilość iteracji
        local it = 0

        # jeśli znaki takie same - nie ma pewności poprawności metody
        if sign(fa) == sign(fb)
            return (r = 0.0, v = 0.0, it = it, err = 1)
        end

        while true
            it = it+1
            # środek przedziału
            local c = (a/2) + (b/2)
            # wartość środka przedziału 
            local fc = f(c)

            # po osiagnieciu precyzji - wyjdź
            if abs(c-a) < delta || abs(fc) < epsilon
                return (r= c, v= fc, it= it, err = 0)
            end

            # wybierz kolejny przedział
            if sign(fc) == sign(fa) 
                fa = fc
                a = c
            else
                fb = fc
                b = c
            end
        end
    end

    #= = =
    # funkcja znajduje zero funkcji metodą stycznych (Newtona)
    #
    # f - funkcja
    # pf - pochodna funkcji
    # x0 - przyblizenie poczatkowe
    # delta - dokładność obliczeń - przedział
    # epsilon - dokładność obliczeń - wartość
    # maxit - maksymalna liczba iteracji
    #
    # return (r: przybliżenie pierwiastka, v: wartość w r, it: ilość iteracji, err: flaga błędu)
    = = =#
    function mstycznych(f, pf, x0::Float64, delta::Float64, epsilon::Float64, maxit:: Int)
        #wartość funkcji w punkcie x0/x1
        local v = f(x0)
        # ilośc iteracji
        local it = 0

        if abs(v) < epsilon
            return (r= x0, v= v, it= it, err = 2)
        end

        while it <= maxit
            it = it+1

            # wartośc pochodnej w x0
            local pf0 = pf(x0)

            # jeśli wartość jest bliska zeru - wyjdź (dostalibyśmy Inf)
            if abs(pf0) < eps(Float64) break; end

            # punkt przybliżany metodą stycznych
            local x1 = x0 - v/pf0
            v = f(x1)

            # po osiagnieciu precyzji - wyjdź
            if abs(x1-x0) < delta || abs(v) < epsilon
                return (r=x1,v=v,it=it,err=0)
            end

            # przygotowanie kolejnej iteracji
            x0 = x1
        end

        return (r=x0,v=v,it=it,err=1)
    end

    #= = =
    # funkcja znajduje zero funkcji metodą siecznych (Eulera)
    #
    # f - funkcja
    # x0 - początek przyblizenia
    # x1 - koniec przyblizenia
    # delta - dokładność obliczeń - przedział
    # epsilon - dokładność obliczeń - wartość
    # maxit - maksymalna liczba iteracji
    #
    # return (r: przybliżenie pierwiastka, v: wartość w r, it: ilość iteracji, err: flaga błędu)
    = = =#
    function msiecznych(f, x0::Float64, x1::Float64, delta::Float64, epsilon::Float64, maxit::Int)
        # wartość funkcji w punkcie x0
        local f0 = f(x0)
        # wartość funkcji w punkcie x1
        local f1 = f(x1)
        # ilość iteracji
        local it = 0

        while it <= maxit
            it = it+1

            # jeśli |fa|>|fb|, zamień wartości
            if abs(f0) > abs(f1)
                x0,x1 = x1,x0
                f0,f1 = f1,f0
            end

            # mianownik współczynnika kierunkowego
            local d = f1-f0
            # jeśli współczynnik kierunkowy jest bliski zeru - wyjdź z błędem
            if abs(d) < eps(Float64)
                break
            end

            # współczynnik kierunkowy prostej
            local s = (x1-x0)/d
            
            #zapamiętaj poprzedni punkt bliższy zeru
            x1 = x0
            f1 = f0

            # przybliżenie kolejnego punktu
            x0 = x0 - f0*s
            f0 = f(x0)

            # po osiagnieciu precyzji - wyjdź
            if abs(x1-x0) < delta || abs(f0) < epsilon
                return (r=x0,v=f0,it=it,err=0)
            end
        end

        return (r=x0,v=f0,it=it,err=1)
    end

end

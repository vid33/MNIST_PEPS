function [ A ] = MPS_GenerateRandom_left_variable( N, Dmax, d )
    Dvec = zeros(1, N-1);
    
    A = cell(N,1);
    
    for kk=1:N-1
        if kk < 15 % prevent matlab from calculating huge powers
            if d^(kk) < Dmax
                Dvec(kk) = d^(kk);
            else 
                Dvec(kk) = Dmax;
            end
        else 
            Dvec(kk) = Dmax;
        end
    end
    
    for kk=1:N-1
        if Dvec(N-kk) <= Dmax && kk < 15
            if Dvec(N-kk) > d^(kk)
                Dvec(N-kk) = d^(kk);
            end
        end
    end

    A{1} = 0.5*(randn(Dvec(1), d)+1i*randn(Dvec(1), d))/sqrt(d);
    A{N} = 0.5*(randn(d, Dvec(N-1))+1i*randn(d, Dvec(N-1)))/sqrt(d);
    for kk=2:(N-1)
        Dav = (1/2)*(Dvec(kk) + Dvec(kk-1)); %this rescaling important, in order to have norm approx 1. Danger of blowup or decay exp in N.
 %       A{kk} = 0.5*(randn(Dvec(kk-1),d,Dvec(kk))+1i*randn(Dvec(kk-1),d,Dvec(kk)))/sqrt(Dav);
        A{kk} = 0.5*(randn(Dvec(kk), d, Dvec(kk-1))+1i*randn(Dvec(kk), d, Dvec(kk-1)))/sqrt(Dav);
    end
end


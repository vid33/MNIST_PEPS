function [ A, l_pos ] = generateRandomPEPS( N1, N2, D, d, digit_no )

    %digit_no for l-index

    A = cell(N1,N2);

    %CORNERS
    A{1,1} = 0.5*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(D);
    
    A{1,N2} = 0.5*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(D);
    A{N1,1} = 0.5*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(D);
    A{N1,N2} = 0.5*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(D);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %SIDES MINUS CORNERS
    for kk=2:(N2-1)
        A{1, kk} = 0.5*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(D);
        A{N1, kk} = 0.5*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(D);
    end
    
    for kk=2:(N1-1)
        A{kk, 1} = 0.5*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(D);
        A{kk, N2} = 0.5*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(D);
    end
    
    %INSIDE
    for kk=2:(N1-1)
        for mm =2:(N2-1)
            A{kk, mm} = 0.5*(randn(D, D, D, D, d)+1i*randn(D, D, D, D, d))/sqrt(D);
        end
    end
    
    %initially l tensor in approx. centre
    l_pos = [round(N1/2), round(N2/2)];
    
    A{l_pos(1), l_pos(2)} = 0.5*(randn(D, D, D, D, d, digit_no)+1i*randn(D, D, D, D, d, digit_no))/sqrt(D);

end


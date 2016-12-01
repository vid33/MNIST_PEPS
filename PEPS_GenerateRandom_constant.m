function [ A, l_pos ] = PEPS_GenerateRandom_constant( N1, N2, D, d, digit_no )

    %digit_no for l-index

    A = cell(N1,N2);

    %CORNERS
    A{1,1} = sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(1/2));
    
    A{1,N2} = sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
    A{N1,1} = sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
    A{N1,N2} = sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %SIDES MINUS CORNERS
    randn_down = sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
    randn_up = sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
    for kk=2:(N2-1)
        A{1, kk} = randn_down;
        A{N1, kk} = randn_up;
    end
    
    for kk=2:(N1-1)
        randn_left = sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
        randn_right = sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
        A{kk, 1} = randn_left;
        A{kk, N2} = randn_right;
    end
    
    %INSIDE
    randn_centre = sqrt(0.5)*(randn(D, D, D, D, d)+1i*randn(D, D, D, D, d))/sqrt(d*D^(4/2));
    for kk=2:(N1-1)
        for mm =2:(N2-1)
            A{kk, mm} = randn_centre;
        end
    end
    
    %initially l tensor in approx. centre
    l_pos = [round(N1/2), round(N2/2)];
    
    A{l_pos(1), l_pos(2)} = sqrt(0.5)*(randn(D, D, D, D, d, digit_no)+1i*randn(D, D, D, D, d, digit_no))/sqrt(d*D^(4/2));

end


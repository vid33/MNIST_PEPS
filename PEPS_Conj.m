function [ PEPS ] = PEPS_Conj( PEPS )

    [N1, N2] = size(PEPS);
    
    for kk=1:N1
        for mm=1:N2
            PEPS{kk,mm} = conj(PEPS{kk,mm});
        end 
    end
    
end


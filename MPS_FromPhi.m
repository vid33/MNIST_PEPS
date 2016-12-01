function [ Phi_lin ] = MPS_FromPhi( Phi )

    [~, SAMPLE_NO]  =size(Phi);

    [N, ~, ~] = size(Phi{1});
    
    Phi_lin = cell(N^2, SAMPLE_NO);
    
    for zz = 1:SAMPLE_NO
        for kk=1:N
            for mm = 1:N
                Phi_lin{N*(kk-1)+ mm, zz} = reshape(Phi{zz}(mm, kk, :), 2, 1);
            end
        
        end

    end
    
    
    [Nlin, ~] = size(Phi_lin);
    
    %Embed into a "zeros" D=2 MPS
    
    for zz = 1:SAMPLE_NO
       
        tmp = zeros(2,2);
        tmp(:, 1) = Phi_lin{1, zz};
        Phi_lin{1, zz} =  tmp;
        
        tmp = zeros(2,2);
        tmp(1, :) = Phi_lin{Nlin, zz};
        Phi_lin{Nlin, zz} =  tmp;
        
        for kk=2:Nlin-1
           
            tmp = zeros(2, 2, 2);
            tmp(1, :, 1) = Phi_lin{kk, zz};
            Phi_lin{kk, zz} =  tmp;
             
        end
        
        
    end
    

end


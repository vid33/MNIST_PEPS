function [ MPS_out ] = PEPS_EnvironmentOfBipartition(PEPS, MPS_init, k_begin, k_end, orientation)

    MPS_out = MPS_init;
    
    if strcmp(orientation, 'left') == true
        for kk=(k_begin+1):k_end
            MPS_out = applyMPOtoMPS(MPS_out, PEPS(:, kk), 'left');
        end
    elseif strcmp(orientation, 'right') == true
        for kk=(k_begin-1):-1:k_end
            MPS_out = applyMPOtoMPS(MPS_out, PEPS(:, kk), 'right'); 
        end
        
    end

end


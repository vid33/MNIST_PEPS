function [ envMPS_left, envMPS_right ] = PEPS_PrepareEnvironments( PEPS, envMPS_left, envMPS_right, gauge_type )

    %The gauge_type variable ('right', 'left', or 'all') determines whether we want to use the envirnments
    %in a scan from the left or from the right, or just calculate all the
    %environments. E.g. 'left' means we calculate right envns from right k=N to k=2
    %and left env at k=1. Minimises memory costs.
    
    % At the moment only works, and is needed for, PEPS = W_Phi_PEPS, but
    % easily generalised

    [SAMPLE_NO, N] = size(envMPS_left);
    
    fprintf('Calculating left and right W_Phi_PEPS environments:\n');
    
    
    if strcmp(gauge_type, 'left') == 1
    
        for image_current=1:SAMPLE_NO
            fprintf('%i ', image_current);
            if mod(image_current,10)==0
                fprintf('\n');
            end
            envMPS_right{image_current, N} = PEPS{image_current}(:,N);
            for kk=N:-1:3
                envMPS_right{image_current, kk-1} = ...
                    PEPS_EnvironmentOfBipartition(PEPS{image_current}, ...
                        envMPS_right{image_current, kk}, kk, kk-1, 'right'); 
            end
            envMPS_right{image_current,1} = []; 

            envMPS_left{image_current,1} = PEPS{image_current}(:,1);
            for kk=2:N
                envMPS_left{image_current,kk} = [];
            end
        end
        fprintf('\n');
    elseif  strcmp(gauge_type, 'right') == 1 %All the left env. MPS, up to rhs boundary, right MPS only at N
        for image_current=1:SAMPLE_NO
            fprintf('%i ', image_current);
            if mod(image_current,10)==0
                fprintf('\n');
            end     
            envMPS_left{image_current,1} = PEPS{image_current}(:,1);
            for kk=1:N-2
                envMPS_left{image_current,kk+1} = ...
                    PEPS_EnvironmentOfBipartition( PEPS{image_current}, ...
                        envMPS_left{image_current,kk}, kk, kk+1, 'left' );
            end            
            envMPS_left{image_current,N} = [];

            envMPS_right{image_current, N} = PEPS{image_current}(:,N);
            for kk=1:N-1
                envMPS_right{image_current,kk} = [];
            end
        end
        
    elseif strcmp(gauge_type, 'all') == 1 %calculate all the environments
        for image_current=1:SAMPLE_NO
            fprintf('%i ', image_current);
            if mod(image_current,10)==0
                fprintf('\n');
            end
            envMPS_right{image_current, N} = PEPS{image_current}(:,N);
            for kk=N:-1:3
                envMPS_right{image_current, kk-1} = ...
                    PEPS_EnvironmentOfBipartition(PEPS{image_current}, ...
                        envMPS_right{image_current, kk}, kk, kk-1, 'right'); 
            end
            envMPS_right{image_current,1} = []; 
            envMPS_left{image_current,1} = PEPS(:,1);
            for kk=1:N-2
                envMPS_left{image_current,kk+1} = ...
                    PEPS_EnvironmentOfBipartition( PEPS, ...
                        envMPS_left{image_current,kk}, kk, kk+1, 'left' );
            end            
            envMPS_left{image_current,N} = [];
            
        end
   end
    
    

end


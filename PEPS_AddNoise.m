function [ PEPS, PEPS_dl, env_mps_left, env_mps_right ] = ...
        PEPS_AddNoise(PEPS, env_mps_left, env_mps_right, lambda, var_type)

    %lambda is a constant determining the magnitude of the random
    % contribution, PEPS_dl == double layer PEPS, 'var_type' can be 'real' or 'cpx'

    if nargin == 4
        var_type = 'real';
    end
    
    N = size(PEPS, 1);
    D = size(PEPS{1,1}, 1);
    d = size(PEPS{1,1}, 3);
    
    if strcmp(var_type, 'real');
        for pos_horiz=1:N
            for pos_vert=1:N
                fprintf('At position (%i, %i)\n', pos_horiz, pos_vert);
                if pos_horiz == 1 &&  pos_vert == 1
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == 1 && pos_vert ==  N
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == 1
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, D, d))/sqrt(d*D^(3/2));
                elseif pos_horiz == N &&  pos_vert == 1
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == N && pos_vert ==  N
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == N
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, D, d))/sqrt(d*D^(3/2));
                elseif pos_vert == 1 %rest of cases we are away from left and right boundaries 
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, D, d))/sqrt(d*D^(3/2));
                elseif pos_vert == N %rest of cases we are away from left and right boundaries 
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, D, d))/sqrt(d*D^(3/2));
                else %away from boundaries 
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*(randn(D, D, D, D, d))/sqrt(d*D^(4/2));
                end

                PEPS_dl = PEPS_Overlap( PEPS, PEPS); %recalculate double layer

                % recalculate norm - NOT OPTIMAL BUT OK FOR 8x8
                if pos_horiz <= 4 %recalculate envn. MPS from left up to pos_horiz=4
                    env_mps_left{4} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,1), 1, 4, 'left' );
                elseif pos_horiz >=5 %recalculate envn. MPS from right up to pos_horiz=5
                    env_mps_right{5} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,N), N, 5, 'right' );
                end

                PEPS_norm = MPS_Overlap(env_mps_left{4}, env_mps_right{5}, 'left', 'right');
                fprintf('PEPS norm before update %d\n', PEPS_norm);

                % reset norm to unity
                PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz}/sqrt(PEPS_norm);

                % need to recalculate only at the point of jump from pos_horiz=4
                % to pos_horiz = 5
                PEPS_dl = PEPS_Overlap( PEPS, PEPS);
                if pos_horiz == 4 && pos_vert == N
                    env_mps_left{4} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,1), 1, 4, 'left' );
                end

            end
        end
    elseif strcmp(var_type, 'cpx')
        for pos_horiz=1:N
            for pos_vert=1:N
                fprintf('At position (%i, %i)\n', pos_horiz, pos_vert);
                if pos_horiz == 1 &&  pos_vert == 1
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == 1 && pos_vert ==  N
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == 1
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
                elseif pos_horiz == N &&  pos_vert == 1
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == N && pos_vert ==  N
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, d)+1i*randn(D, D, d))/sqrt(d*D^(2/2));
                elseif pos_horiz == N
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
                elseif pos_vert == 1 %rest of cases we are away from left and right boundaries 
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
                elseif pos_vert == N %rest of cases we are away from left and right boundaries 
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, D, d)+1i*randn(D, D, D, d))/sqrt(d*D^(3/2));
                else %away from boundaries 
                    PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz} + ...
                                                lambda*sqrt(0.5)*(randn(D, D, D, D, d)+ 1i*randn(D, D, D, D, d))/sqrt(d*D^(4/2));
                end

                PEPS_dl = PEPS_Overlap( PEPS, PEPS); %recalculate double layer

                % recalculate norm - NOT OPTIMAL BUT OK FOR 8x8
                if pos_horiz <= 4 %recalculate envn. MPS from left up to pos_horiz=4
                    env_mps_left{4} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,1), 1, 4, 'left' );
                elseif pos_horiz >=5 %recalculate envn. MPS from right up to pos_horiz=5
                    env_mps_right{5} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,N), N, 5, 'right' );
                end

                PEPS_norm = MPS_Overlap(env_mps_left{4}, env_mps_right{5}, 'left', 'right');
                fprintf('PEPS norm before update %d\n', PEPS_norm);

                % reset norm to unity
                PEPS{pos_vert, pos_horiz} = PEPS{pos_vert, pos_horiz}/sqrt(PEPS_norm);

                % need to recalculate only at the point of jump from pos_horiz=4
                % to pos_horiz = 5
                PEPS_dl = PEPS_Overlap( PEPS, PEPS);
                if pos_horiz == 4 && pos_vert == N
                    env_mps_left{4} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,1), 1, 4, 'left' );
                end

            end
        end
        
    else
        error('''var_type'' should be either ''real'' or ''cpx''.'); 
    end
        
    PEPS_dl = PEPS_Overlap( PEPS, PEPS); %recalculate double layer

    % recalculate dl environments one last time
    env_mps_left{4} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,1), 1, 4, 'left' );
    env_mps_right{5} = PEPS_EnvironmentOfBipartition( PEPS_dl, PEPS_dl(:,N), N, 5, 'right' );
    fprintf('Final norm of noisy W_PEPS is %d \n', MPS_Overlap(env_mps_left{4}, env_mps_right{5}, 'left', 'right') );
end


function [ W_PEPS, W_Phi_PEPS, W_Phi_envMPS_left, W_Phi_envMPS_right, W_conjW_envMPS_left, W_conjW_envMPS_right, ...
    W_Phi_overlap, Cfn_start, Cfn_end] = ...
    PEPS_Sweep_r_to_l(W_PEPS, Phi, W_Phi_PEPS, TARGET, dt, W_Phi_envMPS_left, W_Phi_envMPS_right, W_conjW_envMPS_left, ...
    W_conjW_envMPS_right, W_Phi_overlap, training_labels, NORMALISE)
    
    [SAMPLE_NO, N] = size(W_Phi_envMPS_left);

    Cfn_start = 0;
    Cfn_component = zeros(1, SAMPLE_NO);
    
    for pos_horiz=N:-1:1
        fprintf('$$$$$$ COLUMN %i $$$$$$$$$$ \n', pos_horiz);
        
        if pos_horiz < N && pos_horiz > 1 %calculate right environment as we move from right to left.
            parfor image_current=1:SAMPLE_NO
                current_env_right = W_Phi_envMPS_right(image_current, :); % ensures parfor compatibility
                current_env_right{pos_horiz} = ...
                    PEPS_EnvironmentOfBipartition( W_Phi_PEPS{image_current}, ...
                        current_env_right{pos_horiz+1}, pos_horiz+1, pos_horiz, 'right' );
                W_Phi_envMPS_right(image_current, :) = current_env_right;
            end
            % remove W_Phi_envMPS_right{:, pos_horiz)} from memory - no
            % longer needed
            W_Phi_envMPS_left(:, pos_horiz) = {[]};
        end   
        for pos_vert=N:-1:1
            fprintf('At pos (%i, %i)\n', pos_horiz, pos_vert);
            delta_tensor = 0; % variation of tensor at hand
            if pos_horiz == 1
                parfor kk=1:SAMPLE_NO
                    W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, 1}, W_Phi_envMPS_right{kk, 2}, 'left', 'right');
                    Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
                    var_tmp = conj(Cfn_component(kk))*PEPS_EnvironmentOfSite(W_Phi_envMPS_left{kk,1},...
                        W_Phi_envMPS_right{kk,2}, [pos_vert,0]);
                    delta_tensor = delta_tensor + tensorTimesVec(var_tmp, reshape(Phi{kk}(pos_vert, pos_horiz, :), 2, 1));
                end
            elseif pos_horiz == N
                parfor kk=1:SAMPLE_NO
                    current_env_left = W_Phi_envMPS_left(kk, :);
                    current_env_right = W_Phi_envMPS_right(kk, :);
                    W_Phi_overlap(kk) =  MPS_Overlap(current_env_left{N-1}, current_env_right{N}, 'left', 'right');
                    Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
                    var_tmp = conj(Cfn_component(kk))*PEPS_EnvironmentOfSite(current_env_left{N-1}, ...
                        current_env_right{N}, [0,pos_vert]);
                    delta_tensor = delta_tensor + tensorTimesVec(var_tmp, reshape(Phi{kk}(pos_vert, pos_horiz, :), 2, 1));        
                end
            else % away from left or right boundary
                parfor kk=1:SAMPLE_NO
                    current_env_left = W_Phi_envMPS_left(kk, :);
                    current_env_right = W_Phi_envMPS_right(kk, :);
                    W_Phi_overlap(kk) =  MPS_Overlap(current_env_left{pos_horiz-1}, current_env_right{pos_horiz}, 'left', 'right');
                    Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
                    var_tmp = conj(Cfn_component(kk))*PEPS_EnvironmentOfSite(current_env_left{pos_horiz-1}, ...
                        current_env_right{pos_horiz+1}, [pos_vert, pos_vert], W_Phi_PEPS{kk}(:, pos_horiz));  
                    delta_tensor = delta_tensor + tensorTimesVec(var_tmp, reshape(Phi{kk}(pos_vert, pos_horiz, :), 2, 1));
                end
            end
            
            if pos_horiz+pos_vert < 2*N
                Cfn_old = Cfn;
            end
            Cfn = sum(Cfn_component.*conj(Cfn_component));
            if pos_horiz+pos_vert == 2*N
                Cfn_start = Cfn;
            end
            
            fprintf('Value of cost function is %d.\n', Cfn);
            if pos_horiz + pos_vert < 2*N
                fprintf('Change in cost function is %d\n', Cfn_old - Cfn);
                if Cfn_old - Cfn < 0
                    fprintf('$$$$$$$$$$$$$ WARNING Cfn VALUE INCREASING $$$$$$$$$$$$$$$$\n');
                end
            end

            W_PEPS{pos_vert, pos_horiz} = W_PEPS{pos_vert, pos_horiz} - dt*conj(delta_tensor);
            
            % recalculate PEPS envoronmens, and normalise (?) - THIS IS NOT 100% OPTIMAL, BUT
            % IT'S ALSO NOT VERY EXPENSIVE IN THE CURRENT SETUP TO CALCULATE ALL
            % THE WAY FROM THE BOUNDARIES
            W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS); %recalculate double layer
     
            % calculate norm
            if pos_horiz <= 4 %recalculate W_conjW MPS from left up to pos_horiz=4
                W_conjW_envMPS_left{4} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                    W_conjW_PEPS(:,1), 1, 4, 'left' );
            elseif pos_horiz >=5 %recalculate W_conjW MPS from right up to pos_horiz=5
                    W_conjW_envMPS_right{5} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                                W_conjW_PEPS(:,N), N, 5, 'right' );
            end
            PEPS_norm = MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right');
            fprintf('PEPS norm %d\n', PEPS_norm);
            if NORMALISE == true
                % reset norm to unity
                W_PEPS{pos_vert, pos_horiz} = W_PEPS{pos_vert, pos_horiz}/sqrt(PEPS_norm);

                W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS); %recalculate double layer
                if pos_horiz <= 4 %recalculate W_conjW MPS from left up to pos_horiz=4
                    W_conjW_envMPS_left{4} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                        W_conjW_PEPS(:,1), 1, 4, 'left' );
                elseif pos_horiz >=5 %recalculate W_conjW MPS from right up to pos_horiz=5
                        W_conjW_envMPS_right{5} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                                    W_conjW_PEPS(:,N), N, 5, 'right' );
                end
                PEPS_norm2 = MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right');
                fprintf('PEPS normalisation check %d\n', PEPS_norm2);
            end
            
            % recalculate W_Phi_PEPS, updating only site in question.
            parfor kk=1:SAMPLE_NO
                W_Phi_PEPS{kk}{pos_vert, pos_horiz} = PEPS_OverlapUpdateSite(W_PEPS, Phi{kk}, pos_horiz, pos_vert);
            end
            % and the left environment (or right if we are at pos_horiz=N)
            if pos_horiz == 1
                parfor kk=1:SAMPLE_NO
                    W_Phi_envMPS_left{kk,1} = W_Phi_PEPS{kk}(:,1);
                end
            elseif pos_horiz == N
                parfor kk=1:SAMPLE_NO
                    W_Phi_envMPS_right{kk, N} = W_Phi_PEPS{kk}(:,N);
                end
            else % recalculate from previous position
                parfor kk=1:SAMPLE_NO
                    current_env_right = W_Phi_envMPS_right(kk, :);
                    current_env_right{pos_horiz} = ...
                        PEPS_EnvironmentOfBipartition( W_Phi_PEPS{kk}, ...
                            current_env_right{pos_horiz+1}, pos_horiz+1, pos_horiz, 'right' );
                    W_Phi_envMPS_right(kk, :) = current_env_right;
                end
            end  
        end

    end
    
    Cfn_end = Cfn;
    

end


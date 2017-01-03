function [ W_PEPS, W_Phi_PEPS, W_Phi_envMPS_left, W_Phi_envMPS_right, W_conjW_envMPS_left, W_conjW_envMPS_right, ...
    W_Phi_overlap, Cfn_start, Cfn_end] = ...
    PEPS_Sweep_l_to_r(W_PEPS, Phi, W_Phi_PEPS, TARGET, dt, W_Phi_envMPS_left, W_Phi_envMPS_right, W_conjW_envMPS_left, ...
    W_conjW_envMPS_right, W_Phi_overlap, training_labels)
    
    [SAMPLE_NO, N] = size(W_Phi_envMPS_left);

    Cfn_start = 0;
    Cfn_component = zeros(1, SAMPLE_NO);
    
    for pos_horiz=1:N
        fprintf('$$$$$$ COLUMN %i $$$$$$$$$$ \n', pos_horiz);
        
        if pos_horiz >1 && pos_horiz < N %calculate left environment as we move from left to right.
            for image_current=1:SAMPLE_NO
                W_Phi_envMPS_left{image_current, pos_horiz} = ...
                    PEPS_EnvironmentOfBipartition( W_Phi_PEPS{image_current}, ...
                        W_Phi_envMPS_left{image_current, pos_horiz-1}, pos_horiz-1, pos_horiz, 'left' );
                %remove W_Phi_envMPS_right{:, pos_horiz)} from memory - no
                %longer needed
                W_Phi_envMPS_right{image_current, pos_horiz} = [];
            end
        end   
        for pos_vert=1:N
            fprintf('At pos (%i, %i)\n', pos_horiz, pos_vert);
            delta_tensor = 0; %variation of tensor at hand
            for kk=1:SAMPLE_NO            
                if pos_horiz == 1
                    W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, 1}, W_Phi_envMPS_right{kk, 2}, 'left', 'right');
                    Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
                    var_tmp = conj(Cfn_component(kk))*PEPS_EnvironmentOfSite(W_Phi_envMPS_left{kk,1},...
                        W_Phi_envMPS_right{kk,2}, [pos_vert,0]);
                    delta_tensor = delta_tensor + tensorTimesVec(var_tmp, reshape(Phi{kk}(pos_vert, pos_horiz, :), 2, 1));
                elseif pos_horiz == N
                    W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, N-1}, W_Phi_envMPS_right{kk, N}, 'left', 'right');
                    Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
                    var_tmp = conj(Cfn_component(kk))*PEPS_EnvironmentOfSite(W_Phi_envMPS_left{kk,N-1}, ...
                        W_Phi_envMPS_right{kk,N}, [0,pos_vert]);
                    delta_tensor = delta_tensor + tensorTimesVec(var_tmp, reshape(Phi{kk}(pos_vert, pos_horiz, :), 2, 1));
                else % away from left or right boundary
                    W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, pos_horiz}, W_Phi_envMPS_right{kk, pos_horiz+1}, 'left', 'right');
                    Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
                    var_tmp = conj(Cfn_component(kk))*PEPS_EnvironmentOfSite(W_Phi_envMPS_left{kk,pos_horiz-1}, ...
                        W_Phi_envMPS_right{kk,pos_horiz+1}, [pos_vert, pos_vert], W_Phi_PEPS{kk}(:, pos_horiz));  
                    delta_tensor = delta_tensor + tensorTimesVec(var_tmp, reshape(Phi{kk}(pos_vert, pos_horiz, :), 2, 1));
                end
            end
            
            if pos_horiz+pos_vert > 2
                Cfn_old = Cfn;
            end
            Cfn = sum(Cfn_component.*conj(Cfn_component));
            if pos_horiz+pos_vert == 2
                Cfn_start = Cfn;
            end
            
            fprintf('Value of cost function is %d.\n', Cfn);
            if pos_horiz + pos_vert > 2
                fprintf('Change in cost function is %d\n', Cfn_old - Cfn);
                if Cfn_old - Cfn < 0
                    fprintf('$$$$$$$$$$$$$ WARNING Cfn VALUE INCREASING $$$$$$$$$$$$$$$$\n');
                end
            end
            
            
            % recalculate W_W_ PEPS envoronments, and renormalise - THIS IS NOT 100% OPTIMAL, BUT
            % IT'S ALSO NOT VERY EXPENSIVE for 8x8 TO CALCULATE ALL
            % THE WAY FROM THE BOUNDARIES
            W_PEPS{pos_vert, pos_horiz} = W_PEPS{pos_vert, pos_horiz} - dt*conj(delta_tensor);
            W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS); %recalculate double layer

            % recalculate norm
            if pos_horiz <= 4 %recalculate W_conjW MPS from left up to pos_horiz=4
                W_conjW_envMPS_left{4} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                    W_conjW_PEPS(:,1), 1, 4, 'left' );
            elseif pos_horiz >=5 %recalculate W_conjW MPS from right up to pos_horiz=5
                    W_conjW_envMPS_right{5} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                                W_conjW_PEPS(:,N), N, 5, 'right' );
            end
            PEPS_norm = MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right');
            fprintf('PEPS norm before update %d\n', PEPS_norm);

            %reset norm to unity
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
            fprintf('PEPS normalised check %d\n', PEPS_norm2);

            % recalculate W_Phi_PEPS, updating only site in question.
            for kk=1:SAMPLE_NO
                W_Phi_PEPS{kk}{pos_vert, pos_horiz} = PEPS_OverlapUpdateSite(W_PEPS, Phi{kk}, pos_horiz, pos_vert);
            end
            % and the left environment (or right if we are at pos_horiz=N)
            if pos_horiz == 1
                for kk=1:SAMPLE_NO
                    W_Phi_envMPS_left{kk,1} = W_Phi_PEPS{kk}(:,1);
                end
            elseif pos_horiz == N
                for kk=1:SAMPLE_NO
                    W_Phi_envMPS_right{kk, N} = W_Phi_PEPS{kk}(:,N);
                end
            else %recalculate from previous position
                for kk=1:SAMPLE_NO
                    W_Phi_envMPS_left{kk,pos_horiz} = ...
                        PEPS_EnvironmentOfBipartition( W_Phi_PEPS{kk}, ...
                            W_Phi_envMPS_left{kk, pos_horiz-1}, pos_horiz-1, pos_horiz, 'left' );
                end
            end  
        end

    end
    
    Cfn_end = Cfn;
    

end


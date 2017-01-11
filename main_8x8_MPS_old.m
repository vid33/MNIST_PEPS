clear;

%%%%%  8x8 MPS algorithm   %%%%% 

LOAD_STATE = 'file'; %Possibilites 'file', 'new', or...

D=2; %D MPS
d=2;
N = 8;
N_MPS = N^2;

%The target digit. PEPS will distinguish whether handwriten digit is or is not
%equal to the target value
TARGET = 2;

images_used = 10; %use sum of this many initial images to construct initial PEPS

SAMPLE_NO = 100;

W_Phi_overlap = zeros(1, SAMPLE_NO);
%Cost function components (cost fn is sum (C_component)^2 over training set)
Cfn_component = zeros(1, SAMPLE_NO);  %cost fn is given by Cfn = sum(Cfn_component.^2);
correct_prediction = zeros(1, SAMPLE_NO);

W_Phi_env_left = cell(1, SAMPLE_NO);
W_Phi_env_right = cell(1, SAMPLE_NO);

if strcmp(LOAD_STATE, 'file')
    fIn = sprintf('data/MPSdim=%ix%i_sample_no=%i_D=%i_target=%i_latest',...
        N, N, SAMPLE_NO, D, TARGET);
    load(fIn);   
elseif strcmp(LOAD_STATE, 'new')
   
    fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

    load(fIn);

    fIn_MPS_init = sprintf('data/initial_MPS_N=%i_d=%i_images_used=%i_D_reduced=%i', N, d, images_used, D);
    load(fIn_MPS_init, 'MPS_init');

    W_MPS = MPS_init;
    %W_MPS = MPS_GenerateRandom_left( N, D, d )
    
    W_W_env_left =  MPS_Environments(W_MPS, W_MPS, 'left');
    W_W_env_right =  MPS_Environments(W_MPS, W_MPS, 'right');

    delta_MPS = 2; %constant of variation
    [W_MPS, W_W_env_left, W_W_env_right] = MPS_AddNoise(W_MPS, W_W_env_left, W_W_env_right, delta_MPS, 'cpx');
    
    Phi_small = Phi(:, 1:images_used);
    [ Phi_MPS ] = MPS_FromPhi( Phi );

    for test_image_no = 1:1:images_used;
        fprintf('Overlap of Phi{%i} with truncated MPS is %d\n', ...
         test_image_no, MPS_Overlap(Phi_MPS(:,test_image_no), W_MPS, 'left', 'left'));
    end

    for kk=1:SAMPLE_NO
        fprintf('%i ', kk);
        if mod(kk,10) == 0
            fprintf('\n');
        end
        %W_Phi_overlap(kk) = MPS_Overlap(MPS_init, Phi_MPS(:, kk), 'left', 'left');
        W_Phi_env_left{kk} =  MPS_Environments(W_MPS, Phi_MPS(:, kk), 'left');
        W_Phi_env_right{kk} =  MPS_Environments(W_MPS, Phi_MPS(:, kk), 'right');
        W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{round(N_MPS/2)}, ...
                                W_Phi_env_right{kk}{round(N_MPS/2)+1} }, {[1, 2], [1,2]});

    end
    fprintf('\n');
    
    fOut = sprintf('data/MPSdim=%ix%i_sample_no=%i_D=%i_target=%i_latest',...
        N, N, SAMPLE_NO, D, TARGET);

    save(fOut);
end

%%%%%%%%%%%%%% Learning rate - this will reset whatever is in LOAD_STATE.
dt = 0.001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Initial norm of W_MPS %d \n', MPS_Overlap(W_MPS, W_MPS, 'left', 'left') );
%Contract({ W_W_env_left{round(N_MPS/2)}, W_W_env_right{round(N_MPS/2)+1} }, {[1, 2], [1,2]})

while 1
    % SWEEP l to r
    for pos_mps=1:N_MPS
        fprintf('At position: %i\n', pos_mps);
        
        delta_A=0;
        parfor kk=1:SAMPLE_NO
            if pos_mps == 1
                W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{pos_mps}, W_Phi_env_right{kk}{pos_mps+1} }, ...
                                                {[1, 2], [1,2]});
                Cfn_component(kk) = W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk));
                delta_A = delta_A + conj(Cfn_component(kk))*Contract({W_Phi_env_right{kk}{pos_mps+1}(:, 1), ...
                                Phi_MPS{pos_mps, kk}(1, :)}, {[-1, 1], [1, -2]});
            elseif pos_mps == N_MPS
                W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{pos_mps-1}, W_Phi_env_right{kk}{pos_mps} }, ...
                                                {[1, 2], [1,2]});
                Cfn_component(kk) = W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk));
                delta_A = delta_A + conj(Cfn_component(kk))*Contract({W_Phi_env_left{kk}{pos_mps-1}(:, 1), ...
                                Phi_MPS{pos_mps, kk}(1, :)}, {[-2, 1], [1, -1]});       
            else
                W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{pos_mps-1}, W_Phi_env_right{kk}{pos_mps} }, ...
                                                {[1, 2], [1,2]});
                Cfn_component(kk) = W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk));
                delta_A = delta_A + conj(Cfn_component(kk))*Contract({W_Phi_env_right{kk}{pos_mps+1}(:, 1), ...
                                Phi_MPS{pos_mps, kk}(1, :, 1), W_Phi_env_left{kk}{pos_mps-1}(:, 1)},...
                                {[-1, 1], [1, -2], [-3]});
            end              
        end
        
        if pos_mps > 1
            Cfn_old = Cfn;
        end
        Cfn = sum(Cfn_component.*conj(Cfn_component));
        if pos_mps == 1
            Cfn_start = Cfn;
        elseif pos_mps == N_MPS
            Cfn_end = Cfn;
        end
        
        fprintf('Value of cost function is %d.\n', Cfn);
        if pos_mps > 1
            fprintf('Change in cost function is %d\n', Cfn_old - Cfn);
            if Cfn_old - Cfn < 0
                fprintf('$$$$$$$$$$$$$ WARNING Cfn VALUE INCREASING $$$$$$$$$$$$$$$$\n');
            end
        end
        
        Cfn = sum(Cfn_component.*conj(Cfn_component));
        fprintf('Value of cost function is %d.\n', Cfn);
        
        W_MPS{pos_mps} = W_MPS{pos_mps} - dt*conj(delta_A);
        MPS_norm = MPS_Overlap(W_MPS, W_MPS, 'left', 'left');
       % W_MPS{pos_mps} = W_MPS{pos_mps}/sqrt(MPS_norm);
        fprintf('MPS norm before update %d\n', MPS_norm);
        
        %Update left environment
        for kk=1:SAMPLE_NO
            if pos_mps == 1
                W_Phi_env_left{kk}{pos_mps} = Contract({ W_MPS{pos_mps}, conj(Phi_MPS{pos_mps, kk}) }, {[-1, 1], [-2, 1]});

            elseif pos_mps == N_MPS
                W_Phi_env_right{kk}{pos_mps} = Contract({ W_MPS{pos_mps}, conj(Phi_MPS{pos_mps, kk}) }, {[1, -1], [1, -2]});
            else
                W_Phi_env_left{kk}{pos_mps} = Contract({W_Phi_env_left{kk}{pos_mps-1}, W_MPS{pos_mps}, ...
                            conj(Phi_MPS{pos_mps, kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});
            end   
        end
        
        %break;
    end
    fprintf('Sweep l-r finished, total change in Cfn is: %d\n', Cfn_start - Cfn_end);
    fprintf('Cfn at the end: %d\n', Cfn_end);
        
    if Cfn_start - Cfn_end < 0
        fprintf('Exiting - Cfn no longer decreasing after sweep.\n');
        break;
    end
    
    % SWEEP r to l
    for pos_mps=N_MPS:-1:1
        fprintf('At position: %i\n', pos_mps);
        
        delta_A=0;
        parfor kk=1:SAMPLE_NO
            if pos_mps == 1
                W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{pos_mps}, W_Phi_env_right{kk}{pos_mps+1} }, ...
                                                {[1, 2], [1,2]});
                Cfn_component(kk) = W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk));
                delta_A = delta_A + conj(Cfn_component(kk))*Contract({W_Phi_env_right{kk}{pos_mps+1}(:, 1), ...
                                Phi_MPS{pos_mps, kk}(1, :)}, {[-1, 1], [1, -2]});
            elseif pos_mps == N_MPS
                W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{pos_mps-1}, W_Phi_env_right{kk}{pos_mps} }, ...
                                                {[1, 2], [1,2]});
                Cfn_component(kk) = W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk));
                delta_A = delta_A + conj(Cfn_component(kk))*Contract({W_Phi_env_left{kk}{pos_mps-1}(:, 1), ...
                                Phi_MPS{pos_mps, kk}(1, :)}, {[-2, 1], [1, -1]});       
            else
                W_Phi_overlap(kk) = Contract({ W_Phi_env_left{kk}{pos_mps}, W_Phi_env_right{kk}{pos_mps+1} }, ...
                                                {[1, 2], [1,2]});
                Cfn_component(kk) = W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk));
                delta_A = delta_A + conj(Cfn_component(kk))*Contract({W_Phi_env_right{kk}{pos_mps+1}(:, 1), ...
                                Phi_MPS{pos_mps, kk}(1, :, 1), W_Phi_env_left{kk}{pos_mps-1}(:, 1)},...
                                {[-1, 1], [1, -2], [-3]});
            end
                    
        end
        
        if pos_mps < N_MPS
            Cfn_old = Cfn;
        end
        Cfn = sum(Cfn_component.*conj(Cfn_component));
        if pos_mps == N_MPS
            Cfn_start = Cfn;
        elseif pos_mps == 1
            Cfn_end = Cfn;
        end
        
        fprintf('Value of cost function is %d.\n', Cfn);
        if pos_mps < N_MPS
            fprintf('Change in cost function is %d\n', Cfn_old - Cfn);
            if Cfn_old - Cfn < 0
                fprintf('$$$$$$$$$$$$$ WARNING Cfn VALUE INCREASING $$$$$$$$$$$$$$$$\n');
            end
        end
        
        Cfn = sum(Cfn_component.*conj(Cfn_component));
        fprintf('Value of cost function is %d.\n', Cfn);
        
        W_MPS{pos_mps} = W_MPS{pos_mps} - dt*conj(delta_A);
        MPS_norm = MPS_Overlap(W_MPS, W_MPS, 'left', 'left');
        %W_MPS{pos_mps} = W_MPS{pos_mps}/sqrt(MPS_norm);
        fprintf('MPS norm before update %d\n', MPS_norm);
        
        %Update left environment
        for kk=1:SAMPLE_NO
            if pos_mps == 1
                W_Phi_env_left{kk}{pos_mps} = Contract({ W_MPS{pos_mps}, conj(Phi_MPS{pos_mps, kk}) }, {[-1, 1], [-2, 1]});

            elseif pos_mps == N_MPS
                W_Phi_env_right{kk}{pos_mps} = Contract({ W_MPS{pos_mps}, conj(Phi_MPS{pos_mps, kk}) }, {[1, -1], [1, -2]});
            else
                W_Phi_env_right{kk}{pos_mps} = Contract({W_Phi_env_right{kk}{pos_mps+1}, W_MPS{pos_mps}, ...
                            conj(Phi_MPS{pos_mps, kk})}, {[1, 2], [1, 3, -1], [2, 3, -2]});
            end   
        end
        
        %break;
    end
    fprintf('Sweep r-l finished, total change in Cfn is: %d\n', Cfn_start - Cfn_end);
    fprintf('Cfn at the end: %d\n', Cfn_end);
        
    if Cfn_start - Cfn_end > 0
        fOut = sprintf('data/MPSdim=%ix%i_sample_no=%i_D=%i_target=%i_latest',...
            N, N, SAMPLE_NO, D, TARGET);
        save(fOut);
    elseif Cfn_start - Cfn_end < 0
        fprintf('Exiting - Cfn no longer decreasing after sweep.\n');
       break;
    end
    
    
end




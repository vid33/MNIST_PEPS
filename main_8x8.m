clear;

%%%%%  8x8 brute force algorithm   %%%%% 

LOAD_STATE = 'new'; %Possibilites 'file', 'new', 'file_add_random' or...

digit_no = 10; %0 to 9

D=2; %Bond dimension PEPS
d=2;
N = 8;
SAMPLE_NO = 1000;

%The target digit. PEPS will distinguish whether handwriten digit is or is not
%equal to the target value
TARGET = 2;

if strcmp(LOAD_STATE, 'file')
    fIn = sprintf('data/PEPSdim=%ix%i_sample_no=%i_D=%i_target=%i_latest_from_conv_mps_not_normalised',...
        N, N, SAMPLE_NO, D, TARGET);
    load(fIn);
elseif strcmp(LOAD_STATE, 'file_add_random') %load from file, and add noise to PEPS
    fIn = sprintf('data/PEPSdim=%ix%i_sample_no=%i_D=%i_target=%i_init_basic_from_conv_mps_not_normalised',...
        N, N, SAMPLE_NO, D, TARGET);
    load(fIn);
    
    %Add noise to PEPS, and specifically, directions perpendicular to embedded
    %PEPS are now non-trivial.
    fprintf('Adding noise to PEPS \n');
    delta_PEPS = 0.5; %determines magnitude of random contribution
    [ W_PEPS, W_conjW_PEPS, W_conjW_envMPS_left, W_conjW_envMPS_right ] = ...
            PEPS_AddNoise(W_PEPS, W_conjW_envMPS_left, W_conjW_envMPS_right, delta_PEPS, 'real');
    
    fprintf('Recalculating W_Phi_PEPS structure:\n');
    for kk=1:SAMPLE_NO
        W_Phi_PEPS{kk} = PEPS_Overlap(W_PEPS, Phi{kk});
    end     
    
    %We are going to start the optimisation at the left boundary - 
    %this calculates right envns from right k=N to k=2 and left env only at k=1
    [ W_Phi_envMPS_left, W_Phi_envMPS_right ] = ...
        PEPS_PrepareEnvironments(W_Phi_PEPS, W_Phi_envMPS_left, W_Phi_envMPS_right, 'left');

    %Check if overlap with Phi is still 
    images_test = 20; %test first twenty images
    for image_current = 1:images_test
    fprintf('Ovarlap of Phi{%i} with W, using noisy PEPS, is %d \n', image_current,...
            MPS_Overlap(W_Phi_envMPS_left{image_current,1}, W_Phi_envMPS_right{image_current,2}, 'left', 'right') );
    
    end
    clearvars images_test;
        
    fOut = sprintf('data/PEPSdim=%ix%i_sample_no=%i_D=%i_target=%i_init_basic_random_from_conv_mps_not_normalised',...
        N, N, SAMPLE_NO, D, TARGET);

    save(fOut);
elseif strcmp(LOAD_STATE, 'new')  %Start from scratch
    %Load training set data (training_labels contains the labels).
    fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);
    load(fIn);
    
    W_Phi_overlap = zeros(1, SAMPLE_NO);
    %Cost function components (cost fn is sum (C_component)^2 over training set)
    Cfn_component = zeros(1, SAMPLE_NO);  %cost fn is given by Cfn = sum(Cfn_component.^2);
    prediction = zeros(1, SAMPLE_NO);
    
    N_MPS = N^2;
    W_Phi_PEPS = cell(1, SAMPLE_NO);

    %We will store networks representing contraction of W with conjW, 
    %and of W with Phi, and associated environment MPS.
    W_conjW_envMPS_left = cell(1,N);
    W_conjW_envMPS_right = cell(1,N);
    W_Phi_envMPS_left = cell(SAMPLE_NO,N);
    W_Phi_envMPS_right = cell(SAMPLE_NO,N);

    images_used = 10; %use sum of this many initial images to construct initial PEPS
   % fIn_MPS_init = sprintf('data/initial_MPS_N=%i_d=%i_images_used=%i_D_reduced=%i', N, d, images_used, D);
    fIn_MPS_init = sprintf('data/MPSdim=%ix%i_sample_no=%i_D=%i_target=%i_latest_not_normalised',...
            N, N, SAMPLE_NO, D, TARGET);
    load(fIn_MPS_init, 'MPS_init');

    [ W_PEPS] = embedMPS_IntoPEPS( MPS_init);

    Phi_small = Phi(:, 1:images_used);
    [ Phi_MPS ] = MPS_FromPhi( Phi_small );

    for test_image_no = 1:1:images_used;
        fprintf('Overlap of Phi{%i} with truncated MPS is %d\n', ...
         test_image_no, MPS_Overlap(Phi_MPS(:,test_image_no), MPS_init, 'left', 'left'));
    end

    % Double layer network
    W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS); 

    % Use these to calculate norm, etc.
    W_conjW_envMPS_left{4} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                    W_conjW_PEPS(:,1), 1, 4, 'left' );
    W_conjW_envMPS_right{5} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
                                W_conjW_PEPS(:,N), N, 5, 'right' );
    fprintf('W_PEPS norm is %d \n', ...
        MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right') );
                                                       
    fprintf('Creating W_Phi_PEPS structure:\n');
    for kk=1:SAMPLE_NO
        W_Phi_PEPS{kk} = PEPS_Overlap(W_PEPS, Phi{kk});
    end    
      
    %We are going to start the optimisation at the left boundary - 
    %this calculates right envns from right k=N to k=2 and left env only at k=1
    [ W_Phi_envMPS_left, W_Phi_envMPS_right ] = ...
        PEPS_PrepareEnvironments(W_Phi_PEPS, W_Phi_envMPS_left, W_Phi_envMPS_right, 'left');

    fOut = sprintf('data/PEPSdim=%ix%i_sample_no=%i_D=%i_target=%i_init_from_conv_mps_not_normalised',...
        N, N, SAMPLE_NO, D, TARGET);

    save(fOut);   
end

%%%%%%%%%%%%%% Learning rate - this will reset whatever is in LOAD_STATE.
dt = 0.001;
NORMALISE = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Cost function
for kk=1:SAMPLE_NO
   W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, 1}, W_Phi_envMPS_right{kk, 2}, 'left', 'right');
   %W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, N-1}, W_Phi_envMPS_right{kk, N}, 'left', 'right'); %if starting sweep on rhs
   Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
   if real(W_Phi_overlap(kk)) > 0.5
        prediction(kk) = TARGET;
   else
        prediction(kk) = 0;
   end
end

prediction_target = prediction(training_labels == TARGET);
prediction_non_target = prediction(training_labels ~= TARGET);

fprintf('Percentage of target values correctly identified: %d.\n', ...
        sum(prediction_target == TARGET)/size(prediction_target,2));
fprintf('Percentage of non-target values correctly identified: %d.\n', ...
        sum(prediction_non_target ~= TARGET)/size(prediction_non_target,2));
    
Cfn = sum(Cfn_component.*conj(Cfn_component));
fprintf('Value of cost function is %d.\n', Cfn);

W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS);
% Calculate initial double layer environments, used to calculate norm.
W_conjW_envMPS_left{4} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
    W_conjW_PEPS(:,1), 1, 4, 'left' );
W_conjW_envMPS_right{5} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
            W_conjW_PEPS(:,N), N, 5, 'right' );

fprintf('Initial norm of W_PEPS %d \n', ...
        MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right') );

    
% Main loop, scan the PEPS, apply local grad. desc., 
% until some convergence criterion is satisfied.  

iteration_no = 1;
while 1
    
    fprintf('Entering sweep no %i\n', iteration_no);
    
    % Sweep from left to right
    [ W_PEPS, W_Phi_PEPS, W_Phi_envMPS_left, W_Phi_envMPS_right, W_conjW_envMPS_left, ...
    W_conjW_envMPS_right, W_Phi_overlap, Cfn_start, Cfn_end] = ...
        PEPS_Sweep_l_to_r(W_PEPS, Phi, W_Phi_PEPS, TARGET, dt, W_Phi_envMPS_left, W_Phi_envMPS_right, ... 
                W_conjW_envMPS_left, W_conjW_envMPS_right, W_Phi_overlap, training_labels, NORMALISE);
    W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS);
            
    fprintf('Loop finished, total change in cost function is: %d\n', Cfn_start - Cfn_end);
    fprintf('Cfn at the end: %d\n', Cfn_end);
        
    if Cfn_start - Cfn_end < 0
        fprintf('Exiting - Cfn no longer decreasing after entire scan.\n');
       break;
    end
            
    % Sweep back from right to left
    [ W_PEPS, W_Phi_PEPS, W_Phi_envMPS_left, W_Phi_envMPS_right, W_conjW_envMPS_left, ...
    W_conjW_envMPS_right, W_Phi_overlap, Cfn_start, Cfn_end] = ...
        PEPS_Sweep_r_to_l(W_PEPS, Phi, W_Phi_PEPS, TARGET, dt, W_Phi_envMPS_left, W_Phi_envMPS_right, ... 
                W_conjW_envMPS_left, W_conjW_envMPS_right, W_Phi_overlap, training_labels, NORMALISE);
    W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS);
    
    fprintf('Loop finished, total change in cost function is: %d\n', Cfn_start - Cfn_end);
    fprintf('Cfn at the end: %d\n', Cfn_end);
    
    %Cost function
    for kk=1:SAMPLE_NO
       W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, 1}, W_Phi_envMPS_right{kk, 2}, 'left', 'right');
       %W_Phi_overlap(kk) =  MPS_Overlap(W_Phi_envMPS_left{kk, N-1}, W_Phi_envMPS_right{kk, N}, 'left', 'right'); %if starting sweep on rhs
       Cfn_component(kk) = (W_Phi_overlap(kk) - isequal(TARGET, training_labels(kk))) ;
       if real(W_Phi_overlap(kk)) > 0.5
            prediction(kk) = TARGET;
       else
            prediction(kk) = 0;
       end
    end

    prediction_target = prediction(training_labels == TARGET);
    prediction_non_target = prediction(training_labels ~= TARGET);

    fprintf('Percentage of target values correctly identified: %d.\n', ...
            sum(prediction_target == TARGET)/size(prediction_target,2));
    fprintf('Percentage of non-target values correctly identified: %d.\n', ...
            sum(prediction_non_target ~= TARGET)/size(prediction_non_target,2));

        
    if Cfn_start - Cfn_end > 0
        fOut = sprintf('data/PEPSdim=%ix%i_sample_no=%i_D=%i_target=%i_latest_from_conv_mps_not_normalised',...
            N, N, SAMPLE_NO, D, TARGET);
        save(fOut);
    else
       fprintf('Exiting - Cfn no longer decreasing after entire scan.\n');
       break;
    end
    
    iteration_no = iteration_no+1;
    
end
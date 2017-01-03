clear;

%testing 8x8 functions, no approximations

digit_no = 10; %0 to 9

D=2; %D PEPS
d=2;
N = 8;
SAMPLE_NO = 100;

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

W_conjW_envMPS_left = cell(1,N);
W_conjW_envMPS_right = cell(1,N);
W_conjW_envMPS_up = cell(1,N);
W_conjW_envMPS_down = cell(1,N);

W_Phi_envMPS_left = cell(SAMPLE_NO,N);
W_Phi_envMPS_right = cell(SAMPLE_NO,N);
W_Phi_envMPS_up = cell(SAMPLE_NO, N);
W_Phi_envMPS_down = cell(SAMPLE_NO, N);

MPS_out = cell(SAMPLE_NO,N);

images_used = 10; %use sum of this many initial images to construct initial PEPS

fIn_MPS_init = sprintf('data/initial_MPS_N=%i_d=%i_images_used=%i_D_reduced=%i', N, d, images_used, D);
load(fIn_MPS_init, 'MPS_init');

[ W_PEPS] = embedMPS_IntoPEPS( MPS_init);

Phi_small = Phi(:, 1:images_used);
[ Phi_MPS ] = MPS_FromPhi( Phi_small );

for test_image_no = 1:1:images_used;
    fprintf('Overlap of Phi{%i} with truncated MPS is %d\n', ...
     test_image_no, MPS_Overlap(Phi_MPS(:,test_image_no), MPS_init, 'left', 'left'));
end

%Create networks representing contraction of W with conjW, and of W with
%Phi

W_conjW_PEPS = PEPS_Overlap( W_PEPS, W_PEPS);

%%%%%%% W PEPS norm %%%%%%%%%

W_conjW_envMPS_left{1} = W_conjW_PEPS(:,1);
for kk=1:N-2
    if kk <= 4
        W_conjW_envMPS_left{kk+1} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
            W_conjW_envMPS_left{kk}, kk, kk+1, 'left' );
    else %bond dimension too large
       W_conjW_envMPS_left{kk+1} = []; 
    end
end
W_conjW_envMPS_left{N} = [];

W_conjW_envMPS_right{N} = W_conjW_PEPS(:,N);
for kk=N:-1:3
    if kk >= 5
        W_conjW_envMPS_right{kk-1} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
            W_conjW_envMPS_right{kk}, kk, kk-1, 'right' );
    else %bond dimension too large
       W_conjW_envMPS_right{kk-1} = []; 
    end
end
W_conjW_envMPS_right{1} = [];

fprintf('Norm of W_PEPS is %d \n', ...
        MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right') );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Add noise to PEPS, and specifically, directions perpendicular to embedded
%PEPS are now non-trivial.
delta_PEPS = 0.1; %determines magnitude of random contribution
[ W_PEPS, W_conjW_PEPS, W_conjW_envMPS_left, W_conjW_envMPS_right ] = ...
        PEPS_AddNoise(W_PEPS, W_conjW_envMPS_left, W_conjW_envMPS_right, delta_PEPS);

%%%%%%% W Phi overlap test %%%%%%%%%

for image_current = 1:images_used
    W_Phi_PEPS = PEPS_Overlap(W_PEPS, Phi{image_current});

    W_Phi_envMPS_left{image_current,1} = W_Phi_PEPS(:,1);
    for kk=1:N-2
        W_Phi_envMPS_left{image_current,kk+1} = ...
            PEPS_EnvironmentOfBipartition( W_Phi_PEPS, ...
                W_Phi_envMPS_left{image_current,kk}, kk, kk+1, 'left' );
    end

    W_Phi_envMPS_left{image_current,N} = 0;

    W_Phi_envMPS_right{image_current, N} = W_Phi_PEPS(:,N);
    for kk=N:-1:3
        W_Phi_envMPS_right{image_current, kk-1} = ...
            PEPS_EnvironmentOfBipartition(W_Phi_PEPS, ...
                W_Phi_envMPS_right{image_current, kk}, kk, kk-1, 'right'); 
    end
    W_Phi_envMPS_right{image_current,N} = 0;

    fprintf('Ovarlap of Phi{%i} with W, using PEPS, is %d \n', image_current,...
            MPS_Overlap(W_Phi_envMPS_left{image_current,2}, W_Phi_envMPS_right{image_current,3}, 'left', 'right') );

end
   break; 


%%%%% Norm test 2

W_conjW_envMPS_left{4} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
            W_conjW_PEPS(:,1), 1, 4, 'left' );
        
W_conjW_envMPS_right{5} = PEPS_EnvironmentOfBipartition( W_conjW_PEPS, ...
            W_conjW_PEPS(:,N), N, 5, 'right' );
        
fprintf('Norm of W_PEPS, obtained directly is %d \n', ...
        MPS_Overlap(W_conjW_envMPS_left{4}, W_conjW_envMPS_right{5}, 'left', 'right') );  


size( PEPS_EnvironmentOfSite(W_Phi_envMPS_left{image_current,1}, W_Phi_envMPS_right{image_current,2}, ...
    [3,0]) )

size( PEPS_EnvironmentOfSite(W_Phi_envMPS_left{image_current,1}, W_Phi_envMPS_right{image_current,3},...
    [3,3], W_Phi_PEPS(:, 3)) )

    

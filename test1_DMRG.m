clear;

%testing a bunch subroutines

digit_no = 10; %0 to 9

D=2;
d=2;
N = 14;
SAMPLE_NO = 1000;
Dmax = 10; %max D in envoronment calculations

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

overlapTot = 0;

samplePoints = 1000;


N_MPS = N^2;


Phi_small = Phi(:, 1:Dmax);

[ Phi_MPS ] = MPS_FromPhi( Phi_small );

%Phi_MPS = Phi_MPS(:, 1:15);

[ Phi_sum_MPS ] = MPS_FromPhiSum( Phi_MPS );

%test: approximate Phi_sum_MPS by a lower dimensional MPS B (D = Dmax)
D_final = 2;

Phi_sum_MPS_truncated = MPS_Truncate_SVD(Phi_sum_MPS, D_final);

fprintf('Norm of truncated MPS is %d\n', MPS_Overlap_left_left(Phi_sum_MPS_truncated, Phi_sum_MPS_truncated));
fprintf('Overlap of truncated with original is %d\n', MPS_Overlap_left_left(Phi_sum_MPS_truncated, Phi_sum_MPS));

fprintf('The truncated minus non-truncated norm is %d.\n', MPS_NormDifference_left(Phi_sum_MPS_truncated, Phi_sum_MPS));


B = MPS_GenerateRandom_left( N_MPS, 2, d );
B_norm = MPS_Overlap_left_left(B, B);
B{1} = B{1}/sqrt(B_norm);

%B = Phi_sum_MPS_truncated;


for kk=1:4
    Phi_sum_MPS = MPS_TransformToRightGauge(Phi_sum_MPS);
    B = MPS_ReduceBondDim_left_sweep_up_DMRG( Phi_sum_MPS, B );
    Phi_sum_MPS = MPS_TransformToLeftGauge(Phi_sum_MPS);
    B = MPS_ReduceBondDim_left_sweep_down_DMRG(Phi_sum_MPS, B );
end

[ W_PEPS] = embedMPS_IntoPEPS( Phi_sum_MPS_truncated);

fprintf('Overlap of Phi{4} with truncated MPS is %d\n', MPS_Overlap_left_left(Phi_MPS(:,4), Phi_sum_MPS_truncated));

%chosing a point beyond peps boundaries equivalent to not having an l-index
l_pos = [30, 30];

W_conjW_PEPS = contractW1_W2( W_PEPS, W_PEPS, l_pos); 
W_Phi_PEPS = contractW_Phi(W_PEPS, Phi{4}, l_pos);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_W_MPS_left = W_conjW_PEPS(:,1);
W_W_MPS_right = MPS_Reshape(W_conjW_PEPS(:,N), 'right', 'left'); 

W_W_MPS_left = applyMPOtoMPS_left(W_W_MPS_left, W_conjW_PEPS(:, 2));
W_W_MPS_left = applyMPOtoMPS_left(W_W_MPS_left, W_conjW_PEPS(:, 3));
W_W_MPS_left = applyMPOtoMPS_left(W_W_MPS_left, W_conjW_PEPS(:, 4));

W_W_MPS_right = applyMPOtoMPS_left(W_W_MPS_right, MPO_Reshape(W_conjW_PEPS(:, N-1), 'right', 'left'));
W_W_MPS_right = applyMPOtoMPS_left(W_W_MPS_right, MPO_Reshape(W_conjW_PEPS(:, N-2), 'right', 'left'));
W_W_MPS_right = applyMPOtoMPS_left(W_W_MPS_right, MPO_Reshape(W_conjW_PEPS(:, N-3), 'right', 'left'));
W_W_MPS_right = MPS_Reshape(W_W_MPS_right, 'left', 'right');


double_layer_overlap_test = MPS_Overlap_left_right(W_W_MPS_left, W_W_MPS_right);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_Phi_MPS_left = W_Phi_PEPS(:,1);
W_Phi_MPS_right = MPS_Reshape(W_Phi_PEPS(:,N), 'right', 'left'); 

W_Phi_MPS_left = applyMPOtoMPS_left(W_Phi_MPS_left, W_Phi_PEPS(:, 2));
W_Phi_MPS_left = applyMPOtoMPS_left(W_Phi_MPS_left, W_Phi_PEPS(:, 3));
W_Phi_MPS_left = applyMPOtoMPS_left(W_Phi_MPS_left, W_Phi_PEPS(:, 4));
W_Phi_MPS_left = applyMPOtoMPS_left(W_Phi_MPS_left, W_Phi_PEPS(:, 5));
W_Phi_MPS_left = applyMPOtoMPS_left(W_Phi_MPS_left, W_Phi_PEPS(:, 6));
W_Phi_MPS_left = applyMPOtoMPS_left(W_Phi_MPS_left, W_Phi_PEPS(:, 7));

W_Phi_MPS_right = applyMPOtoMPS_left(W_Phi_MPS_right, MPO_Reshape(W_Phi_PEPS(:, N-1), 'right', 'left'));
W_Phi_MPS_right = applyMPOtoMPS_left(W_Phi_MPS_right, MPO_Reshape(W_Phi_PEPS(:, N-2), 'right', 'left'));
W_Phi_MPS_right = applyMPOtoMPS_left(W_Phi_MPS_right, MPO_Reshape(W_Phi_PEPS(:, N-3), 'right', 'left'));
W_Phi_MPS_right = applyMPOtoMPS_left(W_Phi_MPS_right, MPO_Reshape(W_Phi_PEPS(:, N-4), 'right', 'left'));
W_Phi_MPS_right = applyMPOtoMPS_left(W_Phi_MPS_right, MPO_Reshape(W_Phi_PEPS(:, N-5), 'right', 'left'));
W_Phi_MPS_right = applyMPOtoMPS_left(W_Phi_MPS_right, MPO_Reshape(W_Phi_PEPS(:, N-6), 'right', 'left'));
W_Phi_MPS_right = MPS_Reshape(W_Phi_MPS_right, 'left', 'right');


single_layer_overlap_test = MPS_Overlap_left_right(W_Phi_MPS_left, W_Phi_MPS_right);






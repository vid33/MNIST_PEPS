clear;

%testing a bunch subroutines

digit_no = 10; %0 to 9

D=2;
d=2;
N = 14;
SAMPLE_NO = 1000;
Dmax = 16; %max D in envoronment calculations

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

overlapTot = 0;

samplePoints = 1000;


N_MPS = N^2;


W_MPS1 = MPS_GenerateRandom_left( N_MPS, Dmax, d );

W_MPS2 = MPS_GenerateRandom_left( N_MPS, Dmax, d );


W_MPS1_overlap = MPS_Overlap_left_left(W_MPS1, W_MPS1);
W_MPS1{1} = W_MPS1{1}/sqrt(W_MPS1_overlap);

W_MPS2_overlap = MPS_Overlap_left_left(W_MPS2, W_MPS2);
W_MPS2{1} = W_MPS2{1}/sqrt(W_MPS2_overlap);

[ Wsum ] = MPS_Sum_left_left(W_MPS1, W_MPS2 );

Wsum_overlap = MPS_Overlap_left_left(Wsum, Wsum);

Phi_small = Phi(:, 1:4);

[ Phi_MPS ] = MPS_FromPhi( Phi_small );

%Phi_MPS = Phi_MPS(:, 1:15);

Phi_MPS_overlap = MPS_Overlap_left_left(Phi_MPS(:,3), Phi_MPS(:,3));

[ Phi_sum_MPS ] = MPS_FromPhiSum( Phi_MPS );

%test: approximate f_MPS_left by a lower dimensional MPS B (D = Dmax)
D_final = 15;
[ B ] = MPS_Truncate( W_MPS1, D_final );
B_norm = MPS_Overlap_left_left(B, B);
B{1}=B{1}/sqrt(B_norm);
B_norm_test = MPS_Overlap_left_left(B, B);


B_W_MPS_overlap = MPS_Overlap_left_left(B, W_MPS1);
B_left = MPS_TransformToLeftGauge(B);
   
W_trunc = MPS_Truncate_SVD(W_MPS1, D_final);

fprintf('Norm of truncated is %d\n', MPS_Overlap_left_left(W_trunc, W_trunc));
fprintf('Overlap of truncated with original is %d\n', MPS_Overlap_left_left(W_trunc, W_MPS1));

fprintf('The truncated minus non-truncated norm is %d.\n', MPS_NormDifference_left(W_MPS1, W_trunc));

B = W_trunc;
for kk=1:10
    W_MPS1 = MPS_TransformToRightGauge(W_MPS1);
    B = MPS_ReduceBondDim_left_sweep_up_DMRG( W_MPS1, B );
    W_MPS1 = MPS_TransformToLeftGauge(W_MPS1);
    B = MPS_ReduceBondDim_left_sweep_down_DMRG(W_MPS1, B );
end

break;


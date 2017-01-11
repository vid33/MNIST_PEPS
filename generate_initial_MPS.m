% generate the initial MPS that is used to create the initial PEPS
% bond dimension reduction from sum of Phi states 

clear;

digit_no = 10; %0 to 9

d=2;
N = 8;
SAMPLE_NO = 1000;

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

N_MPS = N^2;

D_MPS_init = 10; %use sum of this many initial images to construct initial PEPS
Phi_small = Phi(:, 1:D_MPS_init);

[ Phi_MPS ] = MPS_FromPhi( Phi_small );

%Phi_MPS = Phi_MPS(:, 1:15);

[ Phi_sum_MPS ] = MPS_FromPhiSum( Phi_MPS );

%test: approximate Phi_sum_MPS by a lower dimensional MPS B 
D_final = 3;

%MPS_init will be the truncation of the sum of Phi states
MPS_init = MPS_Truncate_SVD(Phi_sum_MPS, D_final);

fprintf('Norm of truncated MPS is %d\n', MPS_Overlap(MPS_init, MPS_init, 'left', 'left'));
fprintf('Overlap of truncated with original is %d\n', MPS_Overlap(MPS_init, Phi_sum_MPS, 'left', 'left'));
fprintf('The truncated minus non-truncated norm is %d.\n', MPS_NormDifference_left(MPS_init, Phi_sum_MPS));

for kk=1:50
    Phi_sum_MPS = MPS_TransformToRightGauge(Phi_sum_MPS);
    MPS_init = MPS_ReduceBondDim_left_sweep_up_DMRG( Phi_sum_MPS, MPS_init );
    Phi_sum_MPS = MPS_TransformToLeftGauge(Phi_sum_MPS);
    MPS_init = MPS_ReduceBondDim_left_sweep_down_DMRG(Phi_sum_MPS, MPS_init );
end

%Set all bond dimensions equal to Dmax. 
MPS_init = MPS_ConstantBondDim( MPS_init );

fOut = sprintf('data/initial_MPS_N=%i_d=%i_images_used=%i_D_reduced=%i', N, d, D_MPS_init, D_final);

save(fOut, 'MPS_init');


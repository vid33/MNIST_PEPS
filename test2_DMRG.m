clear;

%testing a bunch subroutines

digit_no = 10; %0 to 9

D=2;
d=2;
N = 14;
SAMPLE_NO = 1000;
Dmax = 8; %max D in envoronment calculations

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

overlapTot = 0;

samplePoints = 1000;

for kk=1:samplePoints

%l_pos 2-component vector indicating the position of l-index
[ W_PEPS, l_pos ] = generateRandomPEPS(N, N, D, d, digit_no); %PEPS tensors

%contraction with image data
%f_PEPS = contractW_Phi(W_PEPS, Phi{1}, l_pos);

%compute double layer PEPS 
%W_conjW_PEPS = contractW1_W2( W_PEPS, W_PEPS, l_pos); 

W_conjW_PEPS = contractW_Phi(W_PEPS, Phi{1}, l_pos);

%test: apply the MPO to left boundar MPS a few times
f_MPS_left = W_conjW_PEPS(:,1);

%MPO_left2 = W_conjW_PEPS(:, 2);
%f_MPS_left = applyMPOtoMPS_left(f_MPS_left, MPO_left2);

%MPO_left3 = W_conjW_PEPS(:, 3);
%f_MPS_left = applyMPOtoMPS_left(f_MPS_left, MPO_left3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%Right environment
    f_MPS_right = reshapeMPS(W_conjW_PEPS(:,N), 'right', 'left'); 
    
%  MPO_right2 = reshapeMPO(W_conjW_PEPS(:, N-1), 'right', 'left');
%f_MPS_right = applyMPOtoMPS_left(f_MPS_right, MPO_right2);

%  MPO_right3 = reshapeMPO(W_conjW_PEPS(:, N-2), 'right', 'left');
%f_MPS_right = applyMPOtoMPS_left(f_MPS_right, MPO_right3);

%    f_MPS_right = reshapeMPS(f_MPS_right, 'left', 'right'); %reshape back


testOverlap = stateOverlap_left_right(f_MPS_left, f_MPS_right);

leftNorm = stateOverlap_left_left(f_MPS_left, f_MPS_left);

break;

overlapTot = overlapTot + testOverlap;

fprintf('Overlap at sample pt. %i/%i is %d\n', kk, samplePoints, testOverlap);

end

fprintf('Average overlap is %d\n', overlapTot/samplePoints);

break;

%test: approximate f_MPS_left by a lower dimensional MPS B (D = Dmax)
[ B ] = generateRandomMPS_left( N, Dmax, d );

%test: do ten sweeps up and down
for kk=1:10
    B = reduceMPSBondDim_left_sweep_up_DMRG( f_MPS_left, B );
    B = reduceMPSBondDim_left_sweep_down_DMRG( f_MPS_left, B );
end

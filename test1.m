clear;

%testing a bunch subroutines

digit_no = 10; %0 to 9

D=2;
d=2;
N = 14;
SAMPLE_NO = 1000;
Dmax = 4; %max D in envoronment calculations

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

%l_pos 2-component vector indicating the position of l-index
[ W_PEPS, l_pos ] = generateRandomPEPS(N, N, D, d, digit_no); %PEPS tensors

%test: the contraction with image data
f_PEPS = contractW_Phi(W_PEPS, Phi{1}, l_pos);

%test: compute double layer PEPS 
W_conjW_PEPS = contractW1_W2( W_PEPS, conjPEPS(W_PEPS), l_pos); 

%test: apply the MPO to left boundar MPS a few times
f_MPS_left = f_PEPS(:,1);

MPO_left2 = f_PEPS(:, 2);

f_MPS_left = applyMPOtoMPS_left(f_MPS_left, MPO_left2);
MPO_left3 = f_PEPS(:, 3);

f_MPS_left = applyMPOtoMPS_left(f_MPS_left, MPO_left2);


%test: approximate f_MPS_left by a lower dimensional MPS B (D = Dmax)
[ B ] = generateRandomMPS_left( N, Dmax, d );

[ rhol_B_B ] = calculate_rhol_left( B, B );
[ rhor_B_B ] = calculate_rhor_left( B, B );
[ rhol_A_B ] = calculate_rhol_left( f_MPS_left, B );
[ rhor_A_B ] = calculate_rhor_left( f_MPS_left, B );

%test: do ten sweeps up and down
for kk=1:10

    B = reduceMPSBondDim_left_sweep_up( f_MPS_left, B, rhor_B_B, rhor_A_B );

    [ rhol_B_B ] = calculate_rhol_left( B, B );
    [ rhol_A_B ] = calculate_rhol_left( f_MPS_left, B );

    B = reduceMPSBondDim_left_sweep_down( f_MPS_left, B, rhol_B_B, rhol_A_B );

    [ rhor_B_B ] = calculate_rhor_left( B, B );
    [ rhor_A_B ] = calculate_rhor_left( f_MPS_left, B );

end






%Test behavior of norms of MPS with randomly generated tensors

clear;

N=64;
D=2;
d=2;

N_test = 1000;
norm_test = zeros(1, N_test);

for kk=1:N_test
    A = MPS_GenerateRandom_left( N, D, d, 'real' );

    norm_test(kk) = MPS_Overlap(A, A, 'left', 'left');
    
    fprintf('Norm of MPS %d \n', norm_test(kk) );
    
end

fprintf('The mean of the norms is %d.\n', mean(norm_test));
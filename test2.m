clear;

digit_no = 10; %0 to 9

D=2;
d=2;
N = 14;
SAMPLE_NO = 1000; %reduced set
Dmax = 4; %max D in envoronment calculations

fIn = sprintf('data/training_data_N=%i_d=%i_samples=%i', N, d, SAMPLE_NO);

load(fIn);

%l_pos 2-component vector indicating the position of l-index
[ W_PEPS, l_pos ] = generateRandomPEPS(N, N, D, d, digit_no); %PEPS tensors

%contraction of PEPS image data object
f_PEPS = contractW_Phi(W_PEPS, Phi{1}, l_pos);


%%%RIGHT ENVIRONMENT MPS
if l_pos(1) < N

    %Reshape everything from right to left and do the work on left, then
    %reshape back at the end
    f_MPS_left = reshapeMPS(f_PEPS(:,N), 'right', 'left'); 

    for kk=N-1:-1:(l_pos(1)+1)
        kk
        MPO_left = reshapeMPO(f_PEPS(:, kk), 'right', 'left');
        f_MPS_left = applyMPOtoMPS_left(f_MPS_left, MPO_left);
   
        fprintf('At row %i. \n', kk);
    
        if D^(N-kk+1) > Dmax;
            %approximate f_MPS_left by a lower dimensional MPS B (D = Dmax)
            if exist('B') == 0
                [ B ] = generateRandomMPS_left( N, Dmax, d ); %better starting point is previous MPS at Dmax
            end

            %do five sweeps up and down%%%%%%%%%%%%
            [ rhol_B_B ] = calculate_rhol_left( B, B );
            [ rhor_B_B ] = calculate_rhor_left( B, B );
            [ rhol_A_B ] = calculate_rhol_left( f_MPS_left, B );
            [ rhor_A_B ] = calculate_rhor_left( f_MPS_left, B );
            for mm=1:5
                B = reduceMPSBondDim_left_sweep_up( f_MPS_left, B, rhor_B_B, rhor_A_B );
                [ rhol_B_B ] = calculate_rhol_left( B, B );
                [ rhol_A_B ] = calculate_rhol_left( f_MPS_left, B );

                B = reduceMPSBondDim_left_sweep_down( f_MPS_left, B, rhol_B_B, rhol_A_B );

                [ rhor_B_B ] = calculate_rhor_left( B, B );
                [ rhor_A_B ] = calculate_rhor_left( f_MPS_left, B );
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            f_MPS_left = B;
            [Dtest, ~] = size(f_MPS_left{1});
            fprintf('Current bond dimension %i\n', Dtest);
            
        end
    end
    f_MPS_right = reshapeMPS(f_MPS_left, 'left', 'right'); %reshape back

end


%%%LEFT ENVIRONMENT MPS
if l_pos(1) >1

    f_MPS_left = f_PEPS(:,1);
    for kk=2:(l_pos(1)-1)
        MPO_left = f_PEPS(:, kk);
        f_MPS_left = applyMPOtoMPS_left(f_MPS_left, MPO_left);
   
        fprintf('At row %i. \n', kk);
    
        if D^kk > Dmax;
            %approximate f_MPS_left by a lower dimensional MPS B (D = Dmax)
            if exist('B') == 0
                [ B ] = generateRandomMPS_left( N, Dmax, d ); %better starting point is previous MPS at Dmax
            end

            %do five sweeps up and down%%%%%%%%%%%%
            [ rhol_B_B ] = calculate_rhol_left( B, B );
            [ rhor_B_B ] = calculate_rhor_left( B, B );
            [ rhol_A_B ] = calculate_rhol_left( f_MPS_left, B );
            [ rhor_A_B ] = calculate_rhor_left( f_MPS_left, B );
            for mm=1:5
                B = reduceMPSBondDim_left_sweep_up( f_MPS_left, B, rhor_B_B, rhor_A_B );
                [ rhol_B_B ] = calculate_rhol_left( B, B );
                [ rhol_A_B ] = calculate_rhol_left( f_MPS_left, B );

                B = reduceMPSBondDim_left_sweep_down( f_MPS_left, B, rhol_B_B, rhol_A_B );

                [ rhor_B_B ] = calculate_rhor_left( B, B );
                [ rhor_A_B ] = calculate_rhor_left( f_MPS_left, B );
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            f_MPS_left = B;
            [Dtest, ~] = size(f_MPS_left{1});
            fprintf('Current bond dimension %i\n', Dtest);
            
        end
    end

end


%%%%%%%%%%%%%%%%% ENVIRONMENT TENSOR FOR l_pos and l_pos+ (0,1)
%%%%%%%%%%%%%%%%%%Contract with the MPO at l_pos, leaving out l_pos and
%%%%%%%%%%%%%%%%% site above it

l_MPO = f_PEPS(:, l_pos(2)); %MPO of columnt in which l_index resides

env_up = Contract({f_MPS_left{N},  l_MPO{N}}, { [1, -3], [ -1, -2, 1]}); 
env_up = Contract({env_up, f_MPS_right{N}}, { [1, -2, -3], [-1, 1]});
env_up = Contract({env_up}, {[-3, -2, -1]});

for kk = 2:(l_pos(2)-1)
    env_up = Contract({env_up, f_MPS_right{N-kk+1}}, {[-1, -2, 1], [1, -4, -3]});
    env_up = Contract({env_up, l_MPO{N-kk+1}}, {[-1, 1, 2, -4], [1, 2, -3, -2]});
    env_up = Contract({env_up, f_MPS_left{N-kk+1}}, {[1, 2, -2, -3], [1, 2, -1]});
end

env_down = Contract({f_MPS_left{1},  l_MPO{1}}, { [-1, 1], [ -2, -3, 1]}); 
env_down = Contract({env_down, f_MPS_right{1}}, { [-1, -2, 1], [-3, 1]});

for kk = 2:(l_pos(2)-1)
    env_up = Contract({env_up, f_MPS_right{N-kk+1}}, {[-1, -2, 1], [-4, 1, -3]});
    env_up = Contract({env_up, l_MPO{N-kk+1}}, {[-1, 1, 2, -4], [-3, 2, 1, -2]});
    env_up = Contract({env_up, f_MPS_left{N-kk+1}}, {[1, 2, -2, -3], [-1, 2, 1]});
end

env_tot = Contract({env_up, f_MPS_left{l_pos(2)+1}}, {[1, -3, -4], [1, -2, -1]}); 
env_tot = Contract({env_tot, f_MPS_right{l_pos(2)+1}}, {[-1, -2, -3, 1], [1, -5, -4]}); 
env_tot = Contract({f_MPS_left{l_pos(2)}, env_tot, f_MPS_right{l_pos(2)}}, ...
                    {[1, -2, -1], [1, -3, -4, -5, 2], [2, -7, -6]});
                
env_tot = Contract({env_tot, env_down}, {[1, -5, -6, -1, -2, -3, 2], [1, -4, 2]}); 














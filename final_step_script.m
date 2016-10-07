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
    env_down = Contract({env_down, f_MPS_right{kk}}, {[-1, -2, 1], [-4, 1, -3]});
    env_down = Contract({env_down, l_MPO{kk}}, {[-1, 1, 2, -4], [-3, 2, 1, -2]});
    env_down = Contract({env_down, f_MPS_left{kk}}, {[1, 2, -2, -3], [-1, 2, 1]});
end

env_tot = Contract({env_up, f_MPS_left{l_pos(2)+1}}, {[1, -3, -4], [1, -2, -1]}); 
env_tot = Contract({env_tot, f_MPS_right{l_pos(2)+1}}, {[-1, -2, -3, 1], [1, -5, -4]}); 
env_tot = Contract({f_MPS_left{l_pos(2)}, env_tot, f_MPS_right{l_pos(2)}}, ...
                    {[1, -2, -1], [1, -3, -4, -5, 2], [2, -7, -6]});

                
env_tot = Contract({env_tot, env_down}, {[1, -5, -6, -1, -2, -3, 2], [1, -4, 2]}); 
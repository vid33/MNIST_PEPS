function [ envs ] = MPS_Environments(A, B, side)

    N = size(A, 1);
    envs = cell(1, N);

    %side refers to 'left' or 'right' environments

    if strcmp(side, 'left')
        for kk=1:N
            if kk==1
               envs{1} = Contract({ A{1}, conj(B{1}) }, {[-1, 1], [-2, 1]}); 
            elseif kk>1 && kk < N
                envs{kk} = Contract({envs{kk-1}, A{kk}, conj(B{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});
            elseif kk == N
                envs{kk} = 0;
            end
        end
    elseif strcmp(side, 'right')
        for kk=N:-1:1
           if kk == N
                envs{kk} = Contract({ A{N}, conj(B{N}) }, {[1, -1], [1, -2]});
           elseif kk<N && kk>1
                envs{kk} = Contract({ envs{kk+1}, A{kk}, conj(B{kk})}, {[1, 2], [1, 3, -1], [2, 3, -2]});
           elseif kk ==1
               envs{kk} = 0;
           end
            
        end
        
    else
        error('last parameter should either be left or right.');
        
    end

end


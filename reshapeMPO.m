function [ Aout ] = reshapeMPO( A, dir1, dir2 )
%reshape an MPO in 'dir1' orientation to 'dir2' orientation
%dir1, dir2 can be 'left', 'right', 'up', 'down'

    if strcmp(dir1, 'down') || strcmp(dir1, 'up')
        [~, N]  =size(A);
    else
        [N, ~] = size(A);
    end
    
    if strcmp(dir2, 'down') || strcmp(dir1, 'up')
        Aout = cell(1, N);
    else
        Aout = cell(N, 1);
    end
    
    
    if strcmp(dir1,'right') && strcmp(dir2,'left') || strcmp(dir1,'left') && strcmp(dir2,'right')
            Aout{1} = Contract({A{N}}, {[-3, -1, -2]});
            Aout{N} = Contract({A{1}}, {[-2, -3, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-3, -4, -1, -2]});
            end
    end
    
    %implement other possibilities...
            

end


function [ Aout ] = MPS_Reshape( A, dir1, dir2 )
%reshape an MPS in 'dir1' orientation to 'dir2' orientation
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
    
    
    if strcmp(dir1,'down') && strcmp(dir2,'left')
            Aout{1} = Contract({A{N}}, {[-1, -2]});
            Aout{N} = Contract({A{1}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-2, -3, -1]});
            end
    elseif strcmp(dir1,'right') && strcmp(dir2,'left')
            Aout{1} = Contract({A{N}}, {[-1, -2]});
            Aout{N} = Contract({A{1}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-3, -1, -2]});
            end  
    elseif strcmp(dir1,'up') && strcmp(dir2,'left')
            Aout = A;      
    elseif strcmp(dir1,'left') && strcmp(dir2,'down')
            Aout{1} = Contract({A{N}}, {[-1, -2]});
            Aout{N} = Contract({A{1}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-3, -1, -2]});
            end
    elseif strcmp(dir1,'left') && strcmp(dir2,'right')
            Aout{1} = Contract({A{N}}, {[-2, -1]});
            Aout{N} = Contract({A{1}}, {[-1, -2]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-2, -3, -1]});
            end
    elseif strcmp(dir1,'up') && strcmp(dir2,'left')
            Aout = A; %no index gymnastics needed here  
    elseif strcmp(dir1, 'up') && strcmp(dir2, 'right')
            Aout{1} = Contract({A{N}}, {[-2, -1]});
            Aout{N} = Contract({A{1}}, {[-1, -2]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-2, -3, -1]});
            end
    elseif strcmp(dir1, 'right') && strcmp(dir2, 'up')
            Aout{1} = Contract({A{N}}, {[-1, -2]});
            Aout{N} = Contract({A{1}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-3, -1, -2]});
            end
    elseif strcmp(dir1, 'down') && strcmp(dir2, 'right')
            Aout{1} = Contract({A{1}}, {[-2, -1]});
            Aout{N} = Contract({A{N}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{kk}}, {[-3, -1, -2]});
            end
    elseif strcmp(dir1, 'right') && strcmp(dir2, 'down')
            Aout{1} = Contract({A{1}}, {[-2, -1]});
            Aout{N} = Contract({A{N}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{kk}}, {[-2, -3, -1]});
            end
    elseif strcmp(dir1, 'down') && strcmp(dir2, 'up')
            Aout{1} = Contract({A{N}}, {[-2, -1]});
            Aout{N} = Contract({A{1}}, {[-1, -2]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-2, -3, -1]});
            end
    elseif strcmp(dir1, 'up') && strcmp(dir2, 'down')
            Aout{1} = Contract({A{N}}, {[-1, -2]});
            Aout{N} = Contract({A{1}}, {[-2, -1]});
            for kk=2:N-1
                Aout{kk} = Contract({A{N-kk+1}}, {[-2, -3, -1]});
            end
    elseif (strcmp(dir1, 'up') || strcmp(dir1, 'down') || strcmp(dir1, 'left')...
              || strcmp(dir1, 'right')) && strcmp(dir1, dir2)
            warning('Input and output orientations are the same.');
            Aout = A;
    else
            error('Specified input and/or output orientations not recognized.');
    end
    


end


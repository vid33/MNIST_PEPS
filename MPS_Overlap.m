function [ overlap ] = MPS_Overlap(A, B, orientation_A, orientation_B)

    [N, ~] = size(A);    
    
    if nargin == 2 %If unspecified, A and B assumed to be in left orientation
        orientation_A = 'left'; orientation_B = 'left';
    end

    if strcmp(orientation_A, 'left') && strcmp(orientation_B, 'left')
    
        overlap = Contract({ A{1}, conj(B{1}) }, {[-1, 1], [-2, 1]});
        for kk=2:N-1      
            overlap = Contract({overlap, A{kk}, conj(B{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});            
        end

        overlap = Contract({overlap, A{N}, conj(B{N})}, {[1, 2], [3, 1], [3, 2]});
    elseif strcmp(orientation_A, 'left') && strcmp(orientation_B, 'right')
        overlap = Contract({ A{1}, B{1} }, {[-1, 1], [-2, 1]});
        for kk=2:N-1
            overlap = Contract({overlap, A{kk},B{kk}}, {[1, 2], [-1, 3, 1], [-2, 2, 3]});            
        end

        overlap = Contract({overlap, A{N}, B{N} }, {[1, 2], [3, 1], [2, 3]});
    else
        error('Other orientations not implemented yet.\n');
    end
    
        
end


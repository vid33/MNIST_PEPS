function [ overlap ] = overlap_A_B_left( A, B )

    %Both A and B are left orientation MPS

    [N, ~] = size(A);

    overlap = Contract({ A{1}, conj(B{1}) }, {[-1, 1], [-2, 1]});
    for kk=2:N-1
            overlap = Contract({overlap, A{kk}, conj(B{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});            
    end

    overlap = Contract({overlap, A{N}, conj(B{N})}, {[1, 2], [3, 1], [3, 2]});
        
end


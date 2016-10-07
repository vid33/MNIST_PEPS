function [ overlap ] = stateOverlap_left_right( A, B )

    %A is a left and B is a right orientatino MPS

    [N, ~] = size(A);

    overlap = Contract({ A{1}, B{1} }, {[-1, 1], [-2, 1]});
    for kk=2:N-1
            overlap = Contract({overlap, A{kk},B{kk}}, {[1, 2], [-1, 3, 1], [-2, 2, 3]});            
    end

    overlap = Contract({overlap, A{N}, B{N} }, {[1, 2], [3, 1], [2, 3]});
        
end


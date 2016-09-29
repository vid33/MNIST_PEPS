function [ rhor ] = calculate_rhor_vert_left( A, B)
    %This is for vertical MPS; find best way to implement horizontal
    %without duplicating code too much?

    [N, ~] = size(A);
    
    rhor = cell(1, N+1);  % convention rhor{N+1} := rhor{0}
    rhor{N} =1;

    rhor{N-1} = Contract({A{N}, conj(B{N})}, {[1, -1], [1, -2]});  
    
    for kk=1:N-2
        rhor{N-kk-1} = Contract({A{N-kk}, rhor{N-kk}}, {[1, -2, -3], [1, -1]});
        rhor{N-kk-1} = Contract({rhor{N-kk-1}, conj(B{N-kk})}, {[1, 2, -1], [1, 2, -2]});     
    end
    
    rhor{N+1} = Contract({A{1}, rhor{1}}, {[1, -1], [1, -2]});
    rhor{N+1} = Contract({rhor{N+1}, conj(B{1})}, {[1, 2], [2, 1]});
    
end


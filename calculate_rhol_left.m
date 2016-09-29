function [ rhol ] = calculate_rhol_left( A, B)

    [N, ~] = size(A);
    
    %MPS environment for A contracted with conj(B)

    rhol = cell(1, N+1);  % rhol{N+1} := rhol{0}
    rhol{N+1} =1;

    rhol{1} = Contract({A{1}, conj(B{1})}, {[-1, 1], [-2, 1]});
    
    for kk=2:N-1
        rhol{kk} = Contract({rhol{kk-1}, A{kk}}, {[1, -3], [-1, -2, 1]});
        rhol{kk} = Contract({rhol{kk}, conj(B{kk})}, {[-1, 1, 2], [-2, 1, 2]});    
    end
    
    rhol{N} = Contract({rhol{N-1}, A{N}}, {[1, -2], [-1, 1]});
    rhol{N} = Contract({rhol{N}, conj(B{N})}, {[1, 2], [1, 2]});

end


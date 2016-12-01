function [ B ] = MPS_Truncate( A, D_final )

    %brute force truncation of an MPS to D_final. Only works if A has fixed
    %bond dimension 

    [N, ~] = size(A);
    [D, d] = size(A{1});
    
    if D_final >= D
        fprintf('Error. Final bond dimension larger than initial.\n');
        return;
    end
        
    B = cell(N, 1);

    B{1} = A{1}(1:D_final, :)*sqrt(d*D^(1/2))/sqrt(d*D_final^(1/2));
    for kk=2:N-1
        B{kk} = A{kk}(1:D_final, :, 1:D_final)*sqrt(d*D)/sqrt(d*D_final);
    end
    B{N} = A{N}(:, 1:D_final)*sqrt(d*D^(1/2))/sqrt(d*D_final^(1/2));
   
end


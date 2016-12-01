function [ B ] = MPS_ReduceBondDim_left_sweep_up( A, B, rhor_B_B, rhor_A_B )

    %pinv is used everywhere, not optimal.
       
    [N, ~] = size(B);
    
    fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
    fprintf('The error at position %i in sweep is %d.\n', 0, normDifferenceBetweenStates_left(A,B));

    B{1} = Contract({ A{1}, rhor_A_B{1} }, { [1, -2], [1, -1]});
    B{1} = Contract({ B{1}, pinv(rhor_B_B{1}) }, { [ 1, -2], [1, -1] });

    fprintf('The error at position %i in sweep is %d.\n', 1, normDifferenceBetweenStates_left(A,B));

    rhol_B_B = Contract({ B{1}, conj(B{1}) }, { [-1, 1], [-2, 1] });
    rhol_A_B = Contract({ A{1}, conj(B{1}) }, { [-1, 1], [-2, 1] });
    for kk=2:N-1
        
        B{kk} = Contract({rhol_A_B, A{kk}}, { [1, -3], [-1, -2, 1]});
        B{kk} = Contract({B{kk}, rhor_A_B{kk}}, { [1, -2, -3], [1, -1] });
        B{kk} = Contract({ pinv(rhol_B_B), B{kk} }, { [1, -3], [-1, -2, 1] } );
        B{kk} = Contract({B{kk}, pinv(rhor_B_B{kk})}, { [1, -2, -3], [1, -1] });
        
        tmp = Contract( {rhol_B_B, B{kk}}, { [1, -3], [-1, -2, 1] });
        rhol_B_B = Contract( {tmp, conj(B{kk})}, { [-1, 1, 2], [-2, 1, 2]});
       
        tmp = Contract( {rhol_A_B, A{kk}}, { [1, -3], [-1, -2, 1] });
        rhol_A_B = Contract( {tmp, conj(B{kk})}, { [-1, 1, 2], [-2, 1, 2]});
        
        fprintf('The error at position %i in sweep is %d.\n', kk, normDifferenceBetweenStates_left(A,B));
        
    end
    
    B{N} = Contract({rhol_A_B, A{N}}, {[1, -2], [-1, 1]});
    B{N} = Contract({pinv(rhol_B_B), B{N}}, {[1, -2], [-1, 1]});  
    
    fprintf('The error at position %i in sweep is %d.\n', N, normDifferenceBetweenStates_left(A,B));

    
end


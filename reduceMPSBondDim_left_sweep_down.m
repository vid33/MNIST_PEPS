function [ B ] = reduceMPSBondDim_left_sweep_down( A, B, rhol_B_B, rhol_A_B )

    %pinv is used everywhere, not optimal.
    
    [N, ~] = size(B);

    fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
    fprintf('The error at position %i in sweep is %d.\n', N+1, normDifferenceBetweenStates_left(A,B));
    
    B{N} = Contract({ A{N}, rhol_A_B{N-1}  }, { [-1, 1], [1, -2]});
        proj_tmp = pinv(rhol_B_B{N-1})*rhol_B_B{N-1};
    B{N} = Contract({ B{N}, proj_tmp*pinv(rhol_B_B{N-1}) }, { [ -1, 1], [1, -2] });
    
    fprintf('The error at position %i in sweep is %d.\n', N, normDifferenceBetweenStates_left(A,B));

    rhor_B_B = Contract({ B{N}, conj(B{N}) }, { [1, -1], [1, -2] });
    rhor_A_B = Contract({ A{N}, conj(B{N}) }, { [1, -1], [1, -2] });
    for kk=(N-1):-1:2
        
        B{kk} = Contract({rhor_A_B, A{kk}}, { [1, -1], [1, -2, -3]});
        B{kk} = Contract({B{kk}, rhol_A_B{kk-1}}, { [-1, -2, 1], [1, -3] });
        
        proj_tmp = rhor_B_B*pinv(rhor_B_B);
        B{kk} = Contract({ pinv(rhor_B_B)*proj_tmp, B{kk} }, { [1, -1], [1, -2, -3]} );
        
        proj_tmp = pinv(rhol_B_B{kk-1})*rhol_B_B{kk-1};
        B{kk} = Contract({B{kk}, proj_tmp*pinv(rhol_B_B{kk-1})}, {  [-1, -2, 1], [1, -3]  });
        
        tmp = Contract( {rhor_B_B, B{kk}}, { [1, -1], [1, -2, -3] });
        rhor_B_B = Contract( {tmp, conj(B{kk})}, { [1, 2, -1], [1, 2, -2]});
       
        tmp = Contract( {rhor_A_B, A{kk}}, { [1, -1], [1, -2, -3] });
        rhor_A_B = Contract( {tmp, conj(B{kk})}, { [1, 2, -1], [1, 2, -2]});
        
        fprintf('The error at position %i in sweep is %d.\n', kk, normDifferenceBetweenStates_left(A,B));
        
    end
    
    B{1} = Contract({rhor_A_B, A{1}}, {[1, -1], [1, -2]});
    B{1} = Contract({pinv(rhor_B_B), B{1}}, {[1, -1], [1, -2]});  
    
    
    fprintf('The error at position %i in sweep is %d.\n', 1, normDifferenceBetweenStates_left(A,B));

    
end


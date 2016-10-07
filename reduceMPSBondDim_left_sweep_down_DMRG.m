function [ B ] = reduceMPSBondDim_left_sweep_down_DMRG( A, B )

    B = transformToLeftGauge(B);
    [ rhol_A_B ] = calculate_rhol_left(A, B);
    
    [N, ~] = size(B);
    [~, d] = size(B{1});


    fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
    fprintf('The error at position %i in sweep is %d.\n', N+1, stateNormDifference_left(A,B));
    
    B{N} = Contract({ A{N}, rhol_A_B{N-1}  }, { [-1, 1], [1, -2]});

    fprintf('The error at position %i in sweep is %d.\n', N, stateNormDifference_left(A,B));

    %QR stuff%%%%%%%%%%%%%%%%%%%
    B{N} = Contract({B{N}}, {[-2, -1]});
    M = reshape(B{N}, size(B{N},1)*size(B{N},2)/d, d); 
    %fix this, unify into single rq function:
    [m, n] = size(M);
    if m > n
        [~, Q] = rq_m_greater_n(M);
    else
        [~, Q] = rq(M);
    end
    
    B{N} = reshape(Q, size(Q,1), d);
    B{N} = Contract({B{N}}, {[-2, -1]});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rhor_A_B = Contract({ A{N}, conj(B{N}) }, { [1, -1], [1, -2] });
    
    for kk=(N-1):-1:2
        
        B{kk} = Contract({rhor_A_B, A{kk}}, { [1, -1], [1, -2, -3]});
        B{kk} = Contract({B{kk}, rhol_A_B{kk-1}}, { [-1, -2, 1], [1, -3] });
  
        fprintf('The error at position %i in sweep is %d.\n', kk, stateNormDifference_left(A,B));

        %QR stuff%%%%%%%%%%%%%%%%%%%
        B{kk} = Contract({B{kk}}, {[-3, -2, -1]});
        M = reshape(B{kk}, size(B{kk},1), d*size(B{kk},3));        
        %fix this, unify into single rq function:
        [m, n] = size(M);
        if m > n
            [~, Q] = rq_m_greater_n(M);
        else
            [~, Q] = rq(M);
        end
        
        if kk==(N-1) 
            index_up = size(B{kk+1},2);
        else
            index_up = size(B{kk+1},3);
        end
        B{kk} = reshape(Q, size(Q,1)*size(Q,2)/(d*index_up), d, index_up);
        B{kk} = Contract({B{kk}}, {[-3, -2, -1]});      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        tmp = Contract( {rhor_A_B, A{kk}}, { [1, -1], [1, -2, -3] });
        rhor_A_B = Contract( {tmp, conj(B{kk})}, { [1, 2, -1], [1, 2, -2]});       
        
    end
    
    B{1} = Contract({rhor_A_B, A{1}}, {[1, -1], [1, -2]});    
    
    fprintf('The error at position %i in sweep is %d.\n', 1, stateNormDifference_left(A,B));

    %QR stuff%%%%%%%%%%%%%%%%%%%
    B{1} = Contract({B{1}}, {[-2, -1]});
    M = reshape(B{1}, 1, size(B{1},2)*d);
    %fix this, unify into single rq function:
    [m, n] = size(M);
    if m > n
        [R, Q] = rq_m_greater_n(M);
    else
        [R, Q] = rq(M);
    end
    
    B{1} = reshape(R*Q, d, size(Q,2)/d); %NB if leftmost R dropped so state normalised.
    B{1} = Contract({B{1}}, {[-2, -1]});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('The error at position %i in sweep is %d.\n', N, stateNormDifference_left(A,B));
    
end


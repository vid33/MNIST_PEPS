function [ B ] = MPS_ReduceBondDim_left_sweep_up_DMRG( A, B )

    %approximate A with B, sweep up

    B = MPS_TransformToRightGauge(B);
    [ rhor_A_B ] = calculate_rhor_left(A, B);
    

    [N, ~] = size(B);
    [~, d] = size(B{1});
    
    fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
    fprintf('The error at position %i in sweep is %d.\n', 0, MPS_NormDifference_left(A,B));

    B{1} = Contract({ A{1}, rhor_A_B{1} }, { [1, -2], [1, -1]});

    
    %QR stuff%%%%%%%%%%%%%%%%%%%
    B{1} = Contract({B{1}}, {[-2, -1]});
    M = reshape(B{1}, d, size(B{1}, 2));
    [Q, ~] = qr(M, 0);
    B{1} = reshape(Q, d, size(Q,2) );
    B{1} = Contract({B{1}}, {[-2, -1]});  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rhol_A_B = Contract({ A{1}, conj(B{1}) }, { [-1, 1], [-2, 1] });   
    
    for kk=2:N-1     
        
        B{kk} = Contract({rhol_A_B, A{kk}}, { [1, -3], [-1, -2, 1]});
        B{kk} = Contract({B{kk}, rhor_A_B{kk}}, { [1, -2, -3], [1, -1] });
        
      % fprintf('The error at position %i in sweep is %d.\n', kk, MPS_NormDifference_left(A,B));
        
        %QR stuff%%%%%%%%%%%%%%%
        B{kk} = Contract({B{kk}}, {[-3, -2, -1]});
        M = reshape(B{kk}, d*size(B{kk},1), size(B{kk},3) );
        [Q, ~] = qr(M, 0);
        B{kk} = reshape(Q, size(B{kk-1}, 1), d, ( size(Q,1)*size(Q,2)/(d*size(B{kk-1},1)) ) );
        B{kk} = Contract({B{kk}}, {[-3, -2, -1]});
        %%%%%%%%%%%%%%%%%
                
        tmp = Contract( {rhol_A_B, A{kk}}, { [1, -3], [-1, -2, 1] });
        rhol_A_B = Contract( {tmp, conj(B{kk})}, { [-1, 1, 2], [-2, 1, 2]});     
        
        %if mod(kk,10) == 0
        %fprintf('The error at position %i in sweep is %d.\n', kk, MPS_NormDifference_left(A,B));
        %end
        
    end
    
    B{N} = Contract({rhol_A_B, A{N}}, {[1, -2], [-1, 1]});
 
    
    %QR stuff%%%%%%%%%%%%%%%%%%%
    B{N} = Contract({B{N}}, {[-2, -1]});
    M = reshape(B{N}, d*size(B{N},1), 1);
    [Q, R] = qr(M, 0);
    B{N} = reshape(Q*sign(R), size(Q,1)/d, d);
    
    %B{N} = reshape(Q*R, size(Q,1)/d, d); %NB if rightmost R dropped state normalised.
    B{N} = Contract({B{N}}, {[-2, -1]});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    fprintf('The error at position %i in sweep is %d.\n', N, MPS_NormDifference_left(A,B));

    fprintf('Norm of B at N is %d \n', MPS_Overlap(B,B, 'left', 'left'));
    
end


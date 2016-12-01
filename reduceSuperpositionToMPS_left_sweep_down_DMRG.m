function [ B ] = reduceSuperpositionToMPS_left_sweep_down_DMRG( Phi_lin, B )

    [SAMPLES, Nlin] = size(Phi_lin);

    B = transformToLeftGauge(B);
    
    rhol_Phi_B = cell(SAMPLES, Nlin+1);
    rhol_Phi_B_full = cell(1, Nlin+1);

    
    [~, d] = size(B{1});
    
    
    %effective left density matrix for Phi overlap w/ conj(B)
    for zz=1:SAMPLES
    
        rhol_Phi_B{zz, Nlin+1} = 1;
    
        rhol_Phi_B{zz, 1} = Contract({Phi_lin{zz, 1}, conj(B{1})}, {[ 1], [-1, 1]});
    
        size(rhol_Phi_B{zz,1})
        for kk=2:Nlin-1
            kk
            rhol_Phi_B{zz, kk} = Contract({rhol_Phi_B{zz, kk-1}, Phi_lin{zz, kk}}, {[-2,1], [-1,1] });
            rhol_Phi_B{zz, kk} = Contract({rhol_Phi_B{zz, kk}, conj(B{kk})}, {[1, 2], [-1, 1, 2]});    
        end
        
        rhol_Phi_B{zz, Nlin} = Contract({rhol_Phi_B{zz, Nlin-1}, Phi_lin{zz, Nlin}}, {[ -2, 1], [-1, 1]});
        rhol_Phi_B{zz, Nlin} = Contract({rhol_Phi_B{zz, Nlin}, conj(B{Nlin})}, {[1, 2], [1, 2]});
        
    end
    
    %for kk=1:Nlin
    %    rhol_Phi_B_full{kk} = rhol_Phi_B{1,kk};
    %    for zz=2:SAMPLES
    %        rhol_Phi_B_full{kk} = rhol_Phi_B_full{kk} + rhol_Phi_B{zz, kk} ;
    %    end
    %end

  %  fprintf('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
  %  fprintf('The error at position %i in sweep is %d.\n', N+1, stateNormDifference_left(A,B));
    
    B{Nlin} = Contract({rhol_Phi_B{1, Nlin-1}, Phi{1, Nlin}}, {[-2], [-1]});
    for zz=2:SAMPLES
        B{Nlin} = B{Nlin} + Contract({rhol_Phi_B{zz, Nlin-1}, Phi{zz, Nlin}}, {[-2], [-1]})
    end
          
    B{Nlin} = Contract({ A{Nlin}, rhol_Phi_B{Nlin-1}  }, { [-1, 1], [1, -2]});

    fprintf('The error at position %i in sweep is %d.\n', N, stateNormDifference_left(A,B));

    %QR stuff%%%%%%%%%%%%%%%%%%%
    B{Nlin} = Contract({B{Nlin}}, {[-2, -1]});
    M = reshape(B{Nlin}, size(B{Nlin},1)*size(B{Nlin},2)/d, d); 
    %fix this, unify into single rq function:
    [m, n] = size(M);
    if m > n
        [~, Q] = rq_m_greater_n(M);
    else
        [~, Q] = rq(M);
    end
    
    B{Nlin} = reshape(Q, size(Q,1), d);
    B{Nlin} = Contract({B{Nlin}}, {[-2, -1]});
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


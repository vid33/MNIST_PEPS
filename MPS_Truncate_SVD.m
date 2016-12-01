function [ A ] = MPS_Truncate_SVD( A, D_final )

    A = MPS_TransformToRightGauge(A);

    [N, ~] = size(A);
    [~, d] = size(A{1});
    
    
   for kk=1:N
       
        if kk==1
            A{kk} = Contract({A{kk}}, {[-2, -1]});
            M = reshape(A{kk}, d, size(A{kk}, 2));
        elseif kk == N
            A{kk} = Contract({A{kk}}, {[-2, -1]});
            M = reshape(A{kk}, d*size(A{kk},1), 1);
        else
            A{kk} = Contract({A{kk}}, {[-3, -2, -1]});
            M = reshape(A{kk}, d*size(A{kk},1), size(A{kk},3) );
        end
            
       % [Q, R] = qr(M, 0);

        [Q, S, V] = svd(M);

%        fprintf('Sum of diag(S) at position %d is %d \n',kk, sum(diag(S.^2)))

        if kk == 1
            if D_final < size(Q,2);
                D_used = D_final;
                Q = Q(1:end, 1:D_used);
                S = S(1:D_used, 1:end);
                S = S./sqrt(sum(diag(S.^2)));
            else
                D_used = size(Q,2);
            end
            A{kk} = reshape(Q, d, D_used );
            A{kk} = Contract({A{kk}}, {[-2, -1]});
            R = S*V';
        elseif kk == N
            Q = Q(1:end, 1);
            R = S(1,1);
            A{kk} = reshape(Q, size(Q,1)/d, d);
            %A{kk} = reshape(Q*R, size(Q,1)/d, d); %NB if rightmost R dropped state normalised.
            A{kk} = Contract({A{kk}}, {[-2, -1]});
        else
            if D_final < size(S,1)
                D_used = D_final;
            else
                D_used = size(S,1);
            end
            
            Q = Q(1:end, 1:D_used);
            S = S(1:D_used, 1:end);
            S = S./sqrt(sum(diag(S.^2)));
            
            R = S*V';
            
            A{kk} = reshape(Q, size(A{kk-1}, 1), d, ( size(Q,1)*size(Q,2)/(d*size(A{kk-1},1)) ) );
            A{kk} = Contract({A{kk}}, {[-3, -2, -1]});
        end

        if kk < N
            if kk+1 == N
                A{kk+1} = Contract({A{kk+1}}, {[-2, -1]});
                A{kk+1} = Contract({R, A{kk+1}}, {[-1, 1], [1, -2]});
                A{kk+1} = Contract({A{kk+1}}, {[-2, -1]});
            else
                A{kk+1} = Contract({A{kk+1}}, {[-3, -2, -1]});
                A{kk+1} = Contract({R, A{kk+1}}, {[-1, 1], [1, -2, -3]});
                A{kk+1} = Contract({A{kk+1}}, {[-3, -2, -1]});
            end

        end
    end

end


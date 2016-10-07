function [ A ] = transformToLeftGauge( A )

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
            
        [Q, R] = qr(M, 0);

        if kk == 1
            A{kk} = reshape(Q, d, size(Q,2) );
            A{kk} = Contract({A{kk}}, {[-2, -1]});
        elseif kk == N
            A{kk} = reshape(Q*R, size(Q,1)/d, d); %NB if rightmost R dropped so state normalised.
            A{kk} = Contract({A{kk}}, {[-2, -1]});
        else
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


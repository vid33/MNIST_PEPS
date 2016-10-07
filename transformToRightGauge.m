function [ A ] = transformToRightGauge( A )

    [N, ~] = size(A);
    [~, d] = size(A{1});

    for kk=N:-1:1
        
        if kk==1
            A{kk} = Contract({A{kk}}, {[-2, -1]});
            M = reshape(A{kk}, 1, size(A{kk},2)*d);
        elseif kk == N
            A{kk} = Contract({A{kk}}, {[-2, -1]});
            M = reshape(A{kk}, size(A{kk},1)*size(A{kk},2)/d, d);
        else
            A{kk} = Contract({A{kk}}, {[-3, -2, -1]});
            M = reshape(A{kk}, size(A{kk},1), d*size(A{kk},3));
        end
            
        %FIX THIS - rq fns downloaded from internet
        [m, n] = size(M);
        if m > n
            [R, Q] = rq_m_greater_n(M);
        else
            [R, Q] = rq(M);
        end
        
        if kk == 1
            A{kk} = reshape(R*Q, d, size(Q,2)/d); %NB if leftmost R dropped so state normalised.
            A{kk} = Contract({A{kk}}, {[-2, -1]});
        elseif kk == N

            A{kk} = reshape(Q, size(Q,1), d); 
            A{kk} = Contract({A{kk}}, {[-2, -1]});
        else
            if kk==(N-1) 
                index_up = size(A{kk+1},2);
            else
                index_up = size(A{kk+1},3);
            end
            A{kk} = reshape(Q, size(Q,1)*size(Q,2)/(d*index_up), d, index_up);
            A{kk} = Contract({A{kk}}, {[-3, -2, -1]});
        end
        
        if kk > 1
            if kk==2
                A{kk-1} = Contract({A{kk-1}}, {[-2, -1]});
                A{kk-1} = Contract({A{kk-1}, R}, {[-1, 1], [1, -2]});
                A{kk-1} = Contract({A{kk-1}}, {[-2, -1]});
            else
                A{kk-1} = Contract({A{kk-1}}, {[-3, -2, -1]});
                A{kk-1} = Contract({A{kk-1}, R}, {[-1, -2, 1], [1, -3]});
                A{kk-1} = Contract({A{kk-1}}, {[-3, -2, -1]});
            end
        end
    end

end


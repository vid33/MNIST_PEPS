function [ out ] = contractW_Phi( W, Phi, l_pos)


    [N, ~] = size(W);
    
    out = cell(N,N);

    %CORNERS
    if l_pos(1) == 1 && l_pos(2) == 1
        out{1,1} = Contract({W{1,1}, reshape(Phi(1, 1, :), 2, 1) }, {[-1, -2, 1, -3], [ 1]});
    else
        out{1,1} = Contract({W{1,1}, reshape(Phi(1, 1, :), 2, 1) }, {[-1, -2, 1], [ 1]});
    end
    
    if l_pos(1) == 1 && l_pos(2) == N
        out{1,N} = Contract({W{1,N}, reshape(Phi(1, N, :), 2, 1) }, {[-1, -2, 1, -3], [ 1]});
    else
        out{1,N} = Contract({W{1,N}, reshape(Phi(1, N, :), 2, 1) }, {[-1, -2, 1], [ 1]});
    end
    
    if l_pos(1) == N && l_pos(2) == 1
        out{N,1} = Contract({W{N,1}, reshape(Phi(N, 1, :), 2, 1) }, {[-1, -2, 1, -3], [ 1]});
    else
        out{N,1} = Contract({W{N,1}, reshape(Phi(N, 1, :), 2, 1) }, {[-1, -2, 1], [ 1]});
    end
    
    if l_pos(1) == N && l_pos(2) == N
        out{N,N} = Contract({W{N,N}, reshape(Phi(N, N, :), 2, 1) }, {[-1, -2, 1, -3], [ 1]});
    else
        out{N,N} = Contract({W{N,N}, reshape(Phi(N, N, :), 2, 1) }, {[-1, -2, 1], [ 1]});
    end
    
    %SIDES MINUS CORNERS
    for kk=2:(N-1)
        if l_pos(1) == 1 && l_pos(2) == kk
            out{1,kk} = Contract({W{1,kk}, reshape(Phi(1, kk, :), 2, 1) }, {[-1, -2, -3, 1, -4], [ 1]});
        else
            out{1,kk} = Contract({W{1,kk}, reshape(Phi(1, kk, :), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        end
        if l_pos(1) == N && l_pos(2) == kk
            out{N,kk} = Contract({W{N,kk}, reshape(Phi(N, kk, :), 2, 1) }, {[-1, -2, -3, 1, -4], [ 1]});
        else
            out{N,kk} = Contract({W{N,kk}, reshape(Phi(N, kk, :), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        end
    end
    
    for kk=2:(N-1)
        if l_pos(1) == kk && l_pos(2) == 1
            out{kk, 1} = Contract({W{kk,1}, reshape(Phi(kk, 1, :), 2, 1) }, {[-1, -2, -3, 1, -4], [ 1]});
        else
            out{kk,1} = Contract({W{kk, 1}, reshape(Phi(kk, 1, :), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        end
        if l_pos(1) == kk && l_pos(2) == N
            out{kk,N} = Contract({W{kk,N}, reshape(Phi(kk, N, :), 2, 1) }, {[-1, -2, -3, 1, -4], [ 1]});
        else
            out{kk,N} = Contract({W{kk,N}, reshape(Phi(kk, N, :), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        end
    end
    
    %INSIDE
    for kk=2:(N-1)
        for mm=2:(N-1)
            if l_pos(1) == kk && l_pos(2) == mm
                out{kk,mm} = Contract({W{kk,mm}, reshape(Phi(kk, mm, :), 2, 1) }, {[-1, -2, -3, -4, 1, -5], [ 1]});
            else
                out{kk,mm} = Contract({W{kk,mm}, reshape(Phi(kk, mm, :), 2, 1) }, {[-1, -2, -3, -4, 1], [ 1]});
            end
        end
    end
    
end


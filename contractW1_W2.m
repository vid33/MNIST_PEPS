function [ out ] = contractW1_W2( W1, W2, l_pos)

    [N, ~] = size(W1); %assume square PEPS
    digit_no = 10;
    
     if l_pos(1) == 1 && l_pos(2) == 1
        [D1, ~, ~, ~] = size(W1{1,1});
        [D2, ~, ~, ~] = size(W2{1,1});
     else
        %[D1, ~, ~] = size(W1{1,1});
        %[D2, ~, ~] = size(W2{1,1});
     end
    
    out = cell(N,N);

    %CORNERS
    if l_pos(1) == 1 && l_pos(2) == 1
        out{1,1} = Contract({W1{1,1}, conj(W2{1,1}) }, {[-1, -3, 1, -5], [-2, -4, 1, -6 ]});
        out{1,1} = reshape(out{1,1}, D1*D2, D1*D2, digit_no, digit_no);
    else
        [D1, D3, ~] = size(W1{1,1}); [D2, D4, ~] = size(W2{1,1});  
        out{1,1} = Contract({W1{1,1}, conj(W2{1,1}) }, {[-1, -3, 1], [-2, -4, 1]});
        out{1,1} = reshape(out{1,1}, D1*D2, D3*D4);
    end
    
    if l_pos(1) == 1 && l_pos(2) == N
        out{1,N} = Contract({W1{1,N}, conj(W2{1,N}) }, {[-1, -3, 1, -5], [ -2, -4, 1, -6]});
        out{1,N} = reshape(out{1,N}, D1*D2, D1*D2, digit_no, digit_no);
    else
        [D1, D3, ~] = size(W1{1,N}); [D2, D4, ~] = size(W2{1,N});  
        out{1,N} = Contract({W1{1,N}, conj(W2{1,N}) }, {[-1, -3, 1], [ -2, -4, 1]});
        out{1,N} = reshape(out{1,N}, D1*D2, D3*D4);    
    end
    
    if l_pos(1) == N && l_pos(2) == 1
        out{N,1} = Contract({W1{N,1}, conj(W2{N,1}) }, {[-1, -3, 1, -5], [ -2, -4, 1, -6]});
        out{N,1} = reshape(out{N,1}, D1*D2, D1*D2, digit_no, digit_no);
    else
        out{N,1} = Contract({W1{N,1}, conj(W2{N,1}) }, {[-1, -3, 1], [ -2, -4, 1]});
        out{N,1} = reshape(out{N,1}, D1*D2, D1*D2);     
    end
    
    if l_pos(1) == N && l_pos(2) == N
        [D1, D3, ~] = size(W1{1,1}); [D2, D4, ~] = size(W2{1,1}); 
        out{N,N} = Contract({W1{N,N}, conj(W2{N,N}) }, {[-1, -3, 1, -5], [ -2, -4, 1, -6]});
        out{N,N} = reshape(out{N,N}, D1*D2, D3*D4, digit_no, digit_no);
    else
        [D1, D3, ~] = size(W1{1,1}); [D2, D4, ~] = size(W2{1,1}); 
        out{N,N} = Contract({W1{N,N}, conj(W2{N,N}) }, {[-1, -3, 1], [ -2, -4, 1]});
        out{N,N} = reshape(out{N,N}, D1*D2, D3*D4); 
    end
    
    %SIDES MINUS CORNERS
    for kk=2:(N-1)
        if l_pos(1) == 1 && l_pos(2) == kk
            out{1,kk} = Contract({W1{1,kk}, conj(W2{1,kk}) }, {[-1, -3, -5, 1, -7], [ -2, -4, -6, 1, -8]});
            out{1,kk} = reshape(out{1,kk}, D1*D2, D1*D2, D1*D2, digit_no, digit_no);
        else
            [D1, D3, D5, ~] = size(W1{1,kk}); [D2, D4, D6, ~] = size(W2{1,kk}); 
            out{1,kk} = Contract({W1{1,kk}, conj(W2{1,kk}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{1,kk} = reshape(out{1,kk}, D1*D2, D3*D4, D5*D6);
        end
        if l_pos(1) == N && l_pos(2) == kk
            out{N,kk} = Contract({W1{N,kk}, conj(W2{N,kk}) }, {[-1, -3, -5, 1, -7], [ -2, -4, -6, 1, -8]});
            out{N,kk} = reshape(out{N,kk}, D1*D2, D1*D2, D1*D2, digit_no, digit_no);
        else
            [D1, D3, D5, ~] = size(W1{N,kk}); [D2, D4, D6, ~] = size(W2{N,kk}); 
            out{N,kk} = Contract({W1{N,kk}, conj(W2{N,kk}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{N,kk} = reshape(out{N,kk}, D1*D2, D3*D4, D5*D6);
        end
    end
    
    for kk=2:(N-1)
        if l_pos(1) == kk && l_pos(2) == 1
            out{kk,1} = Contract({W1{kk,1}, conj(W2{kk,1}) }, {[-1, -3, -5, 1, -7], [ -2, -4, -6, 1, -8]});
            out{kk,1} = reshape(out{kk,1}, D1*D2, D1*D2, D1*D2, digit_no, digit_no);
        else
            [D1, D3, D5, ~] = size(W1{kk,1}); [D2, D4, D6, ~] = size(W2{kk,1}); 
            out{kk,1} = Contract({W1{kk,1}, conj(W2{kk,1}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{kk,1} = reshape(out{kk,1}, D1*D2, D3*D4, D5*D6);
        end
        if l_pos(1) == kk && l_pos(2) == N
            out{kk,N} = Contract({W1{kk,N}, conj(W2{kk,N}) }, {[-1, -3, -5, 1, -7], [ -2, -4, -6, 1, -8]});
            out{kk,N} = reshape(out{kk,N}, D1*D2, D1*D2, D1*D2, digit_no, digit_no);
        else
            [D1, D3, D5, ~] = size(W1{kk,N}); [D2, D4, D6, ~] = size(W2{kk,N}); 
            out{kk,N} = Contract({W1{kk,N}, conj(W2{kk,N}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{kk,N} = reshape(out{kk,N}, D1*D2, D3*D4, D5*D6);
        end
    end
    
    %INSIDE
    for kk=2:(N-1)
        for mm=2:(N-1)
            if l_pos(1) == kk && l_pos(2) == mm
                out{kk,mm} = Contract({W1{kk,mm}, conj(W2{kk,mm}) }, {[-1, -3, -5, -7, 1, -9], [-2, -4, -6, -8, 1, -10]});
                out{kk,mm} = reshape(out{kk,mm}, D1*D2, D1*D2, D1*D2, D1*D2, digit_no, digit_no);
            else
                [D1, D3, D5, D7, ~] = size(W1{kk,mm}); [D2, D4, D6, D8, ~] = size(W2{kk,mm});
                out{kk,mm} = Contract({W1{kk,mm}, conj(W2{kk,mm}) }, {[-1, -3, -5, -7, 1], [-2, -4, -6, -8, 1]});
                out{kk,mm} = reshape(out{kk,mm}, D1*D2, D3*D4, D5*D6, D7*D8);
            end
        end
    end
    
end


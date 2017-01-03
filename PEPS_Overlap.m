function [ out ] = PEPS_Overlap( W1, W2)

% Calculates overlap of two PEPS, i.e. contracts along the physical index,
% and returning the double layer structure (single layer if one of the
% input PEPS is a product state, as is the case for Phi). This is not a completely 
% generic implementation, but only what we need for MNIST - e.g. at the moment
% only W2 is allowed to be a product state. Also, our Phi is stored as a
% tensor, not as a class object

    [N, ~] = size(W1); %assume square PEPS
    out = cell(N,N); 
    
    if  isa(W2, 'cell') == 1 %W2 not a product state
        %CORNERS
        [D1, D3, ~] = size(W1{1,1}); [D2, D4, ~] = size(W2{1,1});  
        out{1,1} = Contract({W1{1,1}, conj(W2{1,1}) }, {[-1, -3, 1], [-2, -4, 1]});
        out{1,1} = reshape(out{1,1}, D1*D2, D3*D4);

        [D1, D3, ~] = size(W1{1,N}); [D2, D4, ~] = size(W2{1,N});  
        out{1,N} = Contract({W1{1,N}, conj(W2{1,N}) }, {[-1, -3, 1], [ -2, -4, 1]});
        out{1,N} = reshape(out{1,N}, D1*D2, D3*D4);    

        [D1, D3, ~] = size(W1{N,1}); [D2, D4, ~] = size(W2{N,1});
        out{N,1} = Contract({W1{N,1}, conj(W2{N,1}) }, {[-1, -3, 1], [ -2, -4, 1]});
        out{N,1} = reshape(out{N,1}, D1*D2, D3*D4);     

        [D1, D3, ~] = size(W1{N,N}); [D2, D4, ~] = size(W2{N,N}); 
        out{N,N} = Contract({W1{N,N}, conj(W2{N,N}) }, {[-1, -3, 1], [ -2, -4, 1]});
        out{N,N} = reshape(out{N,N}, D1*D2, D3*D4); 

        %SIDES MINUS CORNERS
        for kk=2:(N-1)
            [D1, D3, D5, ~] = size(W1{1,kk}); [D2, D4, D6, ~] = size(W2{1,kk}); 
            out{1,kk} = Contract({W1{1,kk}, conj(W2{1,kk}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{1,kk} = reshape(out{1,kk}, D1*D2, D3*D4, D5*D6);

            [D1, D3, D5, ~] = size(W1{N,kk}); [D2, D4, D6, ~] = size(W2{N,kk}); 
            out{N,kk} = Contract({W1{N,kk}, conj(W2{N,kk}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{N,kk} = reshape(out{N,kk}, D1*D2, D3*D4, D5*D6);
        end

        for kk=2:(N-1)
            [D1, D3, D5, ~] = size(W1{kk,1}); [D2, D4, D6, ~] = size(W2{kk,1}); 
            out{kk,1} = Contract({W1{kk,1}, conj(W2{kk,1}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{kk,1} = reshape(out{kk,1}, D1*D2, D3*D4, D5*D6);

            [D1, D3, D5, ~] = size(W1{kk,N}); [D2, D4, D6, ~] = size(W2{kk,N}); 
            out{kk,N} = Contract({W1{kk,N}, conj(W2{kk,N}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out{kk,N} = reshape(out{kk,N}, D1*D2, D3*D4, D5*D6);
        end

        %INSIDE
        for kk=2:(N-1)
            for mm=2:(N-1)
                [D1, D3, D5, D7, ~] = size(W1{kk,mm}); [D2, D4, D6, D8, ~] = size(W2{kk,mm});
                out{kk,mm} = Contract({W1{kk,mm}, conj(W2{kk,mm}) }, {[-1, -3, -5, -7, 1], [-2, -4, -6, -8, 1]});
                out{kk,mm} = reshape(out{kk,mm}, D1*D2, D3*D4, D5*D6, D7*D8);
            end
        end
    
    elseif isa(W2, 'double') == 1 % W2 is a product state; essentially we have W2=Phi in mind 
        %CORNERS
        out{1,1} = Contract({W1{1,1}, reshape(conj(W2(1, 1, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        out{1,N} = Contract({W1{1,N}, reshape(conj(W2(1, N, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        out{N,1} = Contract({W1{N,1}, reshape(conj(W2(N, 1, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        out{N,N} = Contract({W1{N,N}, reshape(conj(W2(N, N, :)), 2, 1) }, {[-1, -2, 1], [ 1]});

        %SIDES MINUS CORNERS
        for kk=2:(N-1)
            out{1,kk} = Contract({W1{1,kk}, reshape(conj(W2(1, kk, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
            out{N,kk} = Contract({W1{N,kk}, reshape(conj(W2(N, kk, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        end

        for kk=2:(N-1)
            out{kk,1} = Contract({W1{kk, 1}, reshape(conj(W2(kk, 1, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
            out{kk,N} = Contract({W1{kk,N}, reshape(conj(W2(kk, N, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        end

        %INSIDE
        for kk=2:(N-1)
            for mm=2:(N-1)
                out{kk,mm} = Contract({W1{kk,mm}, reshape(conj(W2(kk, mm, :)), 2, 1) }, {[-1, -2, -3, -4, 1], [ 1]});
            end
        end
    end
    
end


function [ out ] = PEPS_OverlapUpdateSite( W1, W2, pos_horiz, pos_vert)
% Like PEPS_Overlap, but updates only a single site
% Output is single tensor at site (pos_horiz, pos_vert)

    [N, ~] = size(W1); %assume square PEPS 
    
    if  isa(W2, 'cell') == 1 %W2 not a product state
        %Corners
        if pos_horiz == 1 && pos_vert == 1
            [D1, D3, ~] = size(W1{1,1}); [D2, D4, ~] = size(W2{1,1});  
            out = Contract({W1{1,1}, conj(W2{1,1}) }, {[-1, -3, 1], [-2, -4, 1]});
            out = reshape(out, D1*D2, D3*D4);
        elseif pos_horiz == 1 && pos_vert == N
            [D1, D3, ~] = size(W1{N,1}); [D2, D4, ~] = size(W2{N,1});  
            out = Contract({W1{N,1}, conj(W2{N,1}) }, {[-1, -3, 1], [ -2, -4, 1]});
            out = reshape(out, D1*D2, D3*D4);  
        elseif pos_horiz == N && pos_vert == 1
            [D1, D3, ~] = size(W1{1,N}); [D2, D4, ~] = size(W2{1,N});
            out = Contract({W1{1,N}, conj(W2{1,N}) }, {[-1, -3, 1], [ -2, -4, 1]});
            out = reshape(out, D1*D2, D3*D4);  
        elseif pos_horiz == N && pos_vert == N
            [D1, D3, ~] = size(W1{N,N}); [D2, D4, ~] = size(W2{N,N}); 
            out = Contract({W1{N,N}, conj(W2{N,N}) }, {[-1, -3, 1], [ -2, -4, 1]});
            out = reshape(out, D1*D2, D3*D4); 
        % Sides, but not corners.    
        elseif pos_horiz == 1
            kk = pos_vert;
            [D1, D3, D5, ~] = size(W1{kk,1}); [D2, D4, D6, ~] = size(W2{kk,1}); 
            out = Contract({W1{kk,1}, conj(W2{kk,1}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out = reshape(out, D1*D2, D3*D4, D5*D6);
        elseif pos_horiz == N
            kk = pos_vert;
            [D1, D3, D5, ~] = size(W1{kk,N}); [D2, D4, D6, ~] = size(W2{kk,N}); 
            out = Contract({W1{kk,N}, conj(W2{kk,N}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out = reshape(out, D1*D2, D3*D4, D5*D6);
        elseif pos_vert == 1
            kk = pos_horiz;
            [D1, D3, D5, ~] = size(W1{1,kk}); [D2, D4, D6, ~] = size(W2{1,kk}); 
            out = Contract({W1{1,kk}, conj(W2{1,kk}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out = reshape(out, D1*D2, D3*D4, D5*D6);
        elseif pos_vert == N
            kk = pos_horiz;
            [D1, D3, D5, ~] = size(W1{N,kk}); [D2, D4, D6, ~] = size(W2{N,kk}); 
            out = Contract({W1{N,kk}, conj(W2{N,kk}) }, {[-1, -3, -5, 1], [ -2, -4, -6, 1]});
            out = reshape(out, D1*D2, D3*D4, D5*D6);
        else %Away from edges
            kk = pos_vert; mm = pos_horiz;
            [D1, D3, D5, D7, ~] = size(W1{kk,mm}); [D2, D4, D6, D8, ~] = size(W2{kk,mm});
            out = Contract({W1{kk,mm}, conj(W2{kk,mm}) }, {[-1, -3, -5, -7, 1], [-2, -4, -6, -8, 1]});
            out = reshape(out, D1*D2, D3*D4, D5*D6, D7*D8);
        end
    
    elseif isa(W2, 'double') == 1 % W2 is a product state; essentially we have W2=Phi in mind 
        %CORNERS
        if pos_horiz == 1 && pos_vert == 1
            out = Contract({W1{1,1}, reshape(conj(W2(1, 1, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        elseif pos_horiz == 1 && pos_vert == N
            out = Contract({W1{N,1}, reshape(conj(W2(N, 1, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        elseif pos_horiz == N && pos_vert == 1
            out = Contract({W1{1,N}, reshape(conj(W2(1, N, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        elseif pos_horiz == N && pos_vert == N
            out = Contract({W1{N,N}, reshape(conj(W2(N, N, :)), 2, 1) }, {[-1, -2, 1], [ 1]});
        elseif pos_horiz == 1 %SIDES MINUS CORNERS   
            kk = pos_vert;
            out = Contract({W1{kk,1}, reshape(conj(W2(kk, 1, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        elseif pos_horiz == N
            kk = pos_vert;
            out = Contract({W1{kk, N}, reshape(conj(W2(kk, N, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        elseif pos_vert == 1
            kk = pos_horiz;
            out = Contract({W1{1, kk}, reshape(conj(W2(1, kk, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        elseif pos_vert == N
            kk = pos_horiz;
            out = Contract({W1{N, kk}, reshape(conj(W2(N, kk, :)), 2, 1) }, {[-1, -2, -3, 1], [ 1]});
        else %ON INSIDE
            kk = pos_horiz; mm = pos_vert;
            out = Contract({W1{mm,kk}, reshape(conj(W2(mm, kk, :)), 2, 1) }, {[-1, -2, -3, -4, 1], [ 1]});
        end
    end
    
end


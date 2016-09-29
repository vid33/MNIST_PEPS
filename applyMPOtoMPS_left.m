function [ MPS_out ] = applyMPOtoMPS_left( MPS, MPO)

    %The tensor containing the l-index is always the one being varied, so
    %we don't have to worry about it in environment calculations.
    
    [N, ~] = size(MPS);  
    [D1, d] = size(MPS{1});
    [D2, ~] = size(MPO{1});
    
    MPS_out = cell(N, 1);

    %ENDPOINTS
    MPS_out{1} = Contract({MPS{1}, MPO{1}}, {[-1, 1], [-2, -3, 1]});
    MPS_out{1} = reshape(MPS_out{1}, D1*D2, d);
        
    MPS_out{N} = Contract({MPS{N}, MPO{N}}, {[1, -2], [-1, -3, 1]});
    MPS_out{N} = reshape(MPS_out{N}, d, D1*D2);
        
    %MIDDLE 
    for kk=2:N-1
        MPS_out{kk} = Contract({MPS{kk}, MPO{kk}}, {[-1, 1, -4], [-2, -3, 1, -5]});
        MPS_out{kk} = reshape(MPS_out{kk}, D1*D2, d, D1*D2);
    end
        
end


function [ MPS_out ] = applyMPOtoMPS_left( MPS, MPO)

    %The tensor containing the l-index is always the one being varied, so
    %we don't have to worry about it in environment calculations.
    
    [N, ~] = size(MPS);  
    [~, d] = size(MPS{1});
    
    MPS_out = cell(N, 1);

    %ENDPOINTS
    MPS_out{1} = Contract({MPS{1}, MPO{1}}, {[-1, 1], [-2, -3, 1]});
    MPS_out{1} = reshape(MPS_out{1}, size(MPS{1},1)*size(MPO{1},1), d);
        
    MPS_out{N} = Contract({MPS{N}, MPO{N}}, {[1, -1], [-3, -2, 1]});
    MPS_out{N} = reshape(MPS_out{N}, size(MPS{N},2)*size(MPO{N},2), d  );
    MPS_out{N} = Contract({MPS_out{N} }, { [-2, -1] });
       
    %MIDDLE 
    for kk=2:N-1
        MPS_out{kk} = Contract({MPS{kk}, MPO{kk}}, {[-1, 1, -5], [-2, -3, -4, 1]});
        MPS_out{kk} = reshape(MPS_out{kk}, size(MPS{kk},1)*size(MPO{kk},1), d, ...
                                            size(MPS{kk},3)*size(MPO{kk},3));
    end
        
end


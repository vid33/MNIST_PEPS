function [ tensor_out ] = tensorTimesVec( tensor_in, vec_in  )
% Contract a tensor with up to 4 indices with a vector.

    tensor_size = size(size(tensor_in), 2);
    
    if tensor_size == 1
        tensor_out = Contract({tensor_in, vec_in}, {[-1], [-2]});
    elseif tensor_size == 2
        tensor_out = Contract({tensor_in, vec_in}, {[-1, -2], [-3]});
    elseif tensor_size == 3
        tensor_out = Contract({tensor_in, vec_in}, {[-1, -2, -3], [-4]});
    elseif tensor_size == 4
        tensor_out = Contract({tensor_in, vec_in},{[-1, -2, -3, -4], [-5]});
    else
        error('Input tensor has too many components.');
    end

end


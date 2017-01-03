function [ out ] = MPS_NormDifference_left( A, B )

    out = MPS_Overlap(A, A, 'left', 'left') - MPS_Overlap(A, B, 'left', 'left') - MPS_Overlap(B, A, 'left', 'left') ...
            + MPS_Overlap(B, B, 'left', 'left');

end


function [ out ] = MPS_NormDifference_left( A, B )

    out = MPS_Overlap_left_left(A, A) - MPS_Overlap_left_left(A, B) - MPS_Overlap_left_left(B, A) ...
            + MPS_Overlap_left_left(B, B);

end


function [ out ] = normDifferenceBetweenStates_left( A, B )

    out = overlap_A_B_left(A, A) - overlap_A_B_left(A, B) - overlap_A_B_left(B, A) ...
            + overlap_A_B_left(B, B);


end


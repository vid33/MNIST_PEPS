function [ out ] = stateNormDifference_left( A, B )

    out = stateOverlap_left_left(A, A) - stateOverlap_left_left(A, B) - stateOverlap_left_left(B, A) ...
            + stateOverlap_left_left(B, B);

end


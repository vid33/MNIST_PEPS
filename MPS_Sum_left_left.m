function [ A_and_B ] = MPS_Sum_left_left( A, B )

    %Both A and B are left orientation MPS

    [N, ~] = size(A);
    
    A_and_B = cell(N, 1);
    
    [~, d, ~] = size(A{2});
    
    [D1_A, ~] = size(A{1});
    [D1_B, ~] = size(B{1});
    
    A_and_B{1} = zeros(D1_A + D1_B, d);
    A_and_B{1}(1:D1_A, :) = A{1};
    A_and_B{1}(D1_A+1:end, :) = B{1};
    
    [~, D2_A] = size(A{end});
    [~, D2_B] = size(B{end});

    A_and_B{end} = zeros(d, D2_A + D2_B);   
    A_and_B{end}(:, 1:D2_A) = A{end};
    A_and_B{end}(:, D2_A+1:end) = B{end};
    
    for kk=2:N-1      
        [D1_A, ~, D2_A] = size(A{kk});
        [D1_B, ~, D2_B] = size(B{kk}); 
        A_and_B{kk} = zeros(D1_A + D1_B, d, D2_A + D2_B); 
        A_and_B{kk}(1:D1_A, :, 1:D2_A) = A{kk};
        A_and_B{kk}(D1_A + 1:end, :, D2_A+1:end) = B{kk};
    end
        
end


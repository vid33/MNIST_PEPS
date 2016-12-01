function [ Phi_sum ] = MPS_FromPhiSum( Phi )

    %Both A and B are left orientation MPS

    %D of Phi_sum is equal to the sample size of Phi
    [N, D] = size(Phi);
    
    Phi_sum = cell(N, 1);
    
    [d, ~]  = size(Phi{1});
    
    Phi_sum{1} = zeros(D, d);
    Phi_sum{N} = zeros(d, D);
    
    for kk=1:D
       Phi_sum{1}(kk, :) = (1/sqrt(D))*Phi{1, kk}(:, 1);
       Phi_sum{N}(:, kk) = Phi{N, kk}(1, :);
    end
    
    for zz=2:N-1
       Phi_sum{zz} = zeros(D, d, D); 
       for kk=1:D
          Phi_sum{zz}(kk,:,kk) = Phi{zz, kk}(1,:,1); 
       end
    end
    

end


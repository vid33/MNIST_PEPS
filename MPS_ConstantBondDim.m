function [ A ] = MPS_ConstantBondDim( A )
    % Variable D MPS -> constant D MPS, adding zero directions if necessary
    
    N = size(A, 1);
    d = size(A{1}, 2);
    Dlist = zeros(1, N);
    
    %Find Dmax
    for kk=1:N
        if kk==1
            Dlist(kk) = size(A{kk}, 1);
        elseif kk==N
            Dlist(kk) = size(A{kk}, 2);
        else
            Dlist(kk) = size(A{kk}, 1);
        end
    end
    
    Dmax = max(Dlist);

    for kk=1:N
       
        if kk==1
            [D1, ~]  = size(A{kk});
            if D1 < Dmax
                tmp = zeros(Dmax, d);
                tmp(1:D1, 1:d) = A{kk};
                A{kk} = tmp;
            end
        elseif kk == N
            [~, D1]  = size(A{kk});
            if D1 < Dmax;
                tmp = zeros(d, Dmax);
                tmp(1:d, 1:D1) = A{kk};
                A{kk} = tmp;
            end
        else
            [D1, ~, D2] = size(A{kk});
            if D1 < Dmax || D2 < Dmax
               tmp = zeros(Dmax, d, Dmax);
               tmp(1:D1, 1:d, 1:D2) = A{kk};
               A{kk} = tmp;
               
            end
        end
            
    end

end


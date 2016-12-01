function [ A] = embedMPS_IntoPEPS( MPS)

    [N_MPS, ~] = size(MPS);
    
    N = sqrt(N_MPS);
    
    %This is the implementation for N even -> implement N odd case
    if mod(N,2) == 1
        fprintf('Embedding of MPS into PEPS for N odd not implemented yet.\n');
        return;
    end
        
    A = cell(N,N);

    %CORNERS
    [I1, I3]  = size(MPS{1});
    A{1,1} = zeros(I1, 2, I3);
    A{1,1}(:, 1, :) = MPS{1};

    A{N,1} = Contract({MPS{N}}, {[-1, -3, -2]});
    
    A{N,N} = Contract({MPS{N*(N-1)+1}}, {[-1, -3, -2]});
    
    tmp = Contract({MPS{N*N}}, {[-2, -1]});
    [I1, I3] = size(tmp);
    A{1,N} = zeros(I1, 2, I3);
    A{1,N}(:, 1, :) = tmp;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %sides minus corners, left and right
    for kk=2:(N-1)
       tmp = Contract({MPS{kk}}, {[-1, -3, -2]});
       [I1, I3, I4] = size(tmp);
       A{kk,1} = zeros(I1, 2, I3, I4);
       A{kk,1}(:,1,:,:) = tmp;
       
       tmp = Contract({MPS{N*(N-1)+kk}}, {[-2, -3, -1]});
       [I1, I2, I4] = size(tmp);
       A{N+1-kk,N} = zeros(I1, I2, 2, I4);
       A{N+1-kk,N}(:,:,1,:) = tmp;  
    end
    
    for kk=1:2:(N-2)
        %top 1 - sides minus corners
        tmp = Contract({MPS{N*kk+1}}, {[-1, -3, -2]});
        [I2, I3, I4] = size(tmp);
        A{N,kk+1} = zeros(2, I2, I3, I4);
        A{N,kk+1}(1,:,:,:) = tmp;
        
        %bottom 1 - sides minus corners
        tmp = Contract({MPS{N*kk+N}}, {[-2, -3, -1]});
        [I1, I2, I4] = size(tmp);
        A{1,kk+1} = zeros(I1, I2, 2, I4);
        A{1,kk+1}(:,:,1,:) = tmp;
        
        %inside 1
        for mm=2:N-1
            tmp = Contract({MPS{N*kk+mm}}, {[-2, -3, -1]});
            [I1, I3, I5] = size(tmp);
            A{N+1-mm, kk+1} = zeros(I1, 2, I3, 2, I5);
            A{N+1-mm, kk+1}(:,1,:,1,:) = tmp;
        end
        
    end
    
    for kk=2:2:(N-1)
        %top 2 - sides minus corners
        tmp = Contract({MPS{N*kk+N}}, {[-1, -3, -2]});
        [I1, I2, I4] = size(tmp);
        A{N,kk+1} = zeros(I1, I2, 2, I4);
        A{N,kk+1}(:,:,1,:) = tmp;
        
        %bottom 1 - sides minus corners
        tmp = Contract({MPS{N*kk+1}}, {[-1, -3, -2]});
        [I1, I3, I4] = size(tmp);
        A{1,kk+1} = zeros(I1, 2, I3, I4);
        A{1,kk+1}(:,1,:,:) = tmp;
        
        %inside 1
        for mm=2:N-1
            tmp = Contract({MPS{N*kk+mm}}, {[-1, -3, -2]});
            [I1, I3, I5] = size(tmp);
            A{mm, kk+1} = zeros(I1, 2, I3, 2, I5);
            A{mm, kk+1}(:,1,:,1,:) = tmp;
        end
        
    end

end


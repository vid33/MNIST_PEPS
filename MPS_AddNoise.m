function [ A, env_left, env_right ] = MPS_AddNoise(A, env_left, env_right, lambda, var_type)

    % lambda is a constant determining the magnitude of the random
    % contribution,'var_type' can be 'real' or 'cpx'

    if nargin == 4
        var_type = 'real';
    end
    
    N = size(A, 2);
    D = size(A{2}, 1);
    d = size(A{2}, 2);
    
    if strcmp(var_type, 'real')
    
        for kk=1:N
            if kk==1
                A{1} = A{1} + lambda*(randn(D, d))/sqrt(d*D^(1/2));
                env_left{1} = Contract({ A{1}, conj(A{1}) }, {[-1, 1], [-2, 1]});
                A{1} = A{1}/sqrt(Contract({env_left{1}, env_right{2}}, {[1, 2], [1, 2]}));
                env_left{1} = Contract({ A{1}, conj(A{1}) }, {[-1, 1], [-2, 1]});
            elseif kk == N
                A{N} = A{N} + lambda*(randn(d, D))/sqrt(d*D^(1/2));
                env_right{N}  = Contract({ A{N}, conj(A{N}) }, {[1, -1], [1, -2]});
                A{N} = A{N}/sqrt(Contract({env_left{N-1}, env_right{N}}, {[1, 2], [1, 2]}));
                env_right =  MPS_Environments(A, A, 'right');
            else
                A{kk} = A{kk} + lamdba*(randn(D,d,D))/sqrt(d*D^(2/2));
                env_left{kk} = Contract({env_left{kk-1}, A{kk}, conj(A{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});
                A{kk} = A{kk}/sqrt(Contract({env_left{kk}, env_right{kk+1}}, {[1, 2], [1, 2]}));
                env_left{kk} = Contract({env_left{kk-1}, A{kk}, conj(A{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});
            end
        end
    
    elseif strcmp(var_type, 'cpx')
        for kk=1:N
            if kk==1
                A{1} = A{1} + lambda*sqrt(0.5)*(randn(D, d)+ 1i*randn(D, d))/sqrt(d*D^(1/2));
                env_left{1} = Contract({ A{1}, conj(A{1}) }, {[-1, 1], [-2, 1]});
                A{1} = A{1}/sqrt(Contract({env_left{1}, env_right{2}}, {[1, 2], [1, 2]}));
                env_left{1} = Contract({ A{1}, conj(A{1}) }, {[-1, 1], [-2, 1]});
            elseif kk == N
                A{N} = A{N} + lambda*sqrt(0.5)*(randn(d, D) + 1i*randn(d, D))/sqrt(d*D^(1/2));
                env_right{N}  = Contract({ A{N}, conj(A{N}) }, {[1, -1], [1, -2]});
                A{N} = A{N}/sqrt(Contract({env_left{N-1}, env_right{N}}, {[1, 2], [1, 2]}));
                env_right =  MPS_Environments(A, A, 'right');
            else
                A{kk} = A{kk} + lamdba*sqrt(0.5)*(randn(D,d,D)+1i*randn(D,d,D))/sqrt(d*D^(2/2));
                env_left{kk} = Contract({env_left{kk-1}, A{kk}, conj(A{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});
                A{kk} = A{kk}/sqrt(Contract({env_left{kk}, env_right{kk+1}}, {[1, 2], [1, 2]}));
                env_left{kk} = Contract({env_left{kk-1}, A{kk}, conj(A{kk})}, {[1, 2], [-1, 3, 1], [-2, 3, 2]});
            end
        end    
    else
        error('''var_type'' should be either ''real'' or ''cpx''.');
    end
end


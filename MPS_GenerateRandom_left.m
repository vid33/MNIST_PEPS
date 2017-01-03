function [ A ] = MPS_GenerateRandom_left( N, D, d, type )
    %type can be real ('real') or complex ('cpx')
    
    if nargin == 3
        type = 'cpx';
    end
    
    A = cell(N,1);

    if strcmp(type, 'cpx')
        A{1} = sqrt(0.5)*(randn(D, d)+1i*randn(D, d))/sqrt(d*D^(1/2));
        A{N} = sqrt(0.5)*(randn(d, D)+1i*randn(d, D))/sqrt(d*D^(1/2));
        for kk=2:(N-1)
            A{kk} = sqrt(0.5)*(randn(D,d,D)+1i*randn(D,d,D))/sqrt(d*D^(2/2));
        end
    elseif strcmp(type, 'real')
        A{1} = (randn(D, d))/sqrt(d*D^(1/2));
        A{N} = (randn(d, D))/sqrt(d*D^(1/2));
        for kk=2:(N-1)
            A{kk} = (randn(D,d,D))/sqrt(d*D^(2/2));
        end  
    else
        error('''type'' parameter must be ''real'' or ''cpx''.');
    end
    
end


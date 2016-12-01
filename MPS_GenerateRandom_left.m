function [ A ] = MPS_GenerateRandom_left( N, D, d )

    A = cell(N,1);

    A{1} = sqrt(0.5)*(randn(D, d)+1i*randn(D, d))/sqrt(d*D^(1/2));
    A{N} = sqrt(0.5)*(randn(d, D)+1i*randn(d, D))/sqrt(d*D^(1/2));
    for kk=2:(N-1)
        A{kk} = sqrt(0.5)*(randn(D,d,D)+1i*randn(D,d,D))/sqrt(d*D^(2/2));
    end
end


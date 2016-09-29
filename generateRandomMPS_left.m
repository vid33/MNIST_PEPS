function [ A ] = generateRandomMPS_left( N, D, d )

    A = cell(N,1);

    A{1} = 0.5*(randn(D, d)+1i*randn(D, d))/sqrt(D);
    A{N} = 0.5*(randn(d, D)+1i*randn(d, D))/sqrt(D);
    for kk=2:(N-1)
        A{kk} = 0.5*(randn(D,d,D)+1i*randn(D,d,D))/sqrt(D);
    end
end


clear;

digit_no = 10;
D=4; d=2;

A1 = 0.5*(randn(D, d, D, digit_no)+1i*randn(D, d, D, digit_no))/sqrt(D);

A2 = 0.5*(randn(D, d, D)+1i*randn(D, d, D))/sqrt(D);


B1 = Contract({A1, A2}, {[-1, -2, 1, -5], [1, -3, -4]});

B2 = reshape(B1, D, d^2, D, digit_no);

B1mat = reshape(B1, d*D, d*D*digit_no);

[U, S, V] = svd(B1mat);

%NB contracted index is now 
C1 = reshape(U, D, d, D*d);

C2 = reshape(S*V', D*d, d, D, digit_no);



B1_test = Contract({C1, C2}, {[-1, -2, 1], [1, -3, -4, -5]});
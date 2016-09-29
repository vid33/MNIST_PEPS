function [ M_pinv ] = pinv_eig( M )
   %Eigenvalue decomp pseudoinverse. 
   
   [D, ~] = size(M);
   
   [V, Diag,] = eig(M);
   
   for kk=1:D
      if abs(Diag(kk,kk)) > 1e-15
          fprintf('Was here\n')
          Diag(kk,kk) = 1/(Diag(kk,kk));
      else
          Diag(kk, kk) = 0;
      end
   end
   
   M_pinv = V*Diag/V;
      
end


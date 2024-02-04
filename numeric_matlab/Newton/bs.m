% Code 12: Backward Substitution for Upper Triangular Systems 
% Input U: Upp. Triangular non-singular square matrix 
% b: column right-hand side 
% Output x: solution of Ux = b 
function x = bs(U,b) 
x = 0*b; n = length(b); x(n) = b(n)/U(n,n); 
for ii = n-1:-1:1 
 x(ii) = (b(ii)-U(ii,ii+1:n)*x(ii+1:n))/U(ii,ii); 
end 
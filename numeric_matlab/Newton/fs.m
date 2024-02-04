% Code 11: Forward Substitution for Lower Triangular Systems 
% Input L: Low Triangular non-singular square matrix 
% b: column right-hand side 
% Output x: solution of Lx = b 
function x = fs(L,b) 
x = 0*b; n = length(b); x(1) = b(1)/L(1,1); 
for ii = 2:n 
 x(ii) = (b(ii)-L(ii,1:ii-1)*x(1:ii-1))/L(ii,ii); 
end 
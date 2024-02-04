% Code 14: PA = LU (Solver for Ax = b) 
% Input: L (unit lower triangular matrix) 
% U (upper triangular matrix) 
% P (reordering vector) 
% b (right-hand side) 
% Output: solution x 
% Note: requires fs.m and bs.m (Codes 11 & 12) 
function x = plusolve(L,U,P,b) 
n = length(b); 
for k = 1:n-1; b([k P(k)]) = b([P(k) k]); end 
y = fs(L,b); x = bs(U,y); 
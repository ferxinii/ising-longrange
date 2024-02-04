% Code 19: Computation of the Jacobian J
% Input: F(x):R^m ---> R^n
% x: (m x 1)-vector ; F: (n x 1)-vector
% Output: DF(x) (n x m) Jacobian matrix at x
function DF = jac(F,x)
    f1 = feval(F,x); n = length(f1); m = length(x);
    DF = zeros(n,m); H = sqrt(eps)*eye(m);
    for j = 1:m
        f2 = feval(F,x+H(:,j)); DF(:,j) = (f2 - f1)/H(j,j);
    end
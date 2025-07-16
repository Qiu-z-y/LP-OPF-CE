% The code can derive the voltage amplitude and phase angle from matrix A. 
% If an error occurs, the truncation error for calculating the null space of matrix A can be appropriately adjusted.

balance_bus = 1;
tol = 1e-3;
[U, S, V] = svd(value(A));
s = diag(S);
null_space_custom = V(:, s < tol);
 X1 = null_space_custom(1:ns/2,1);
 X2 = null_space_custom(ns+1 : ns / 2 * 3,1);
 X3 = null_space_custom(ns / 2  + 1 : ns,1);
 X4 = null_space_custom(ns / 2 * 3 + 1 : 2*ns,1);

syms x y;
eq1 = -X1(balance_bus) * x + X2(balance_bus) * y == 1;
eq2 = X2(balance_bus) * x + X1(balance_bus) * y == 0;
sol1 = solve([eq1, eq2], [x, y]);

syms z w;
eq1 = -X3(balance_bus) * z + X4(balance_bus) * w == 1;
eq2 = X4(balance_bus) * z + X3(balance_bus) * w == 0;
sol2 = solve([eq1, eq2], [z, w]);

X = [ -(X1 + X2*1i).*(double(sol1.x) + double(sol1.y)*1i);  -(X3 + X4*1i).*(double(sol2.z) + double(sol2.w)*1i)];

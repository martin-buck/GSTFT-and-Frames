% Build a given matrix
a = 1/6;
b = -2/3;
G = [1 a a b b a b a a a; ...
    a 1 a b a a a a b b; ...
    a a 1 a a b b a b a; ...
    b b a 1 a b a a a a; ...
    b a a a 1 a a b b a; ...
    a a b b a 1 a b a a; ...
    b a b a a a 1 a a b; ...
    a a a a b b a 1 a b; ...
    a b b a b a a a 1 a; ...
    a b a a a a b b a 1];
% To show it is a Gram matrix for a frame consisting of N vectors in R^d:
% 1) G^2=AG, where A is a constant (not a matrix)
% 2) G^*=G
% 3) G is positive semi-definite

% 1) Check if the above polynomial matrix equation is satisfied
G_squared = G*G;
[V, D] = eig(G);
lambda = D(end);
disp(G_squared-lambda*G)

% 2) Check if G is symmetric
sym = issymmetric(G);
disp(sym)

% 3) Check if the eigenvalues are not negative
disp(D)


A = [0 1 0 0 1; 1 0 1 0 0; 0 1 0 1 0; 0 0 1 0 1; 1 0 0 1 0];
D = diag([2 2 2 2 2]);
L = D - A;

[V1, D1] = eig(L);

[V2, D2] = eig(A);

Lswrt = sqrtm(L);

dft = dftmtx(5);
idft = dft';

% compute the rho_G
rho = 4;

d1_sqrt = sqrt(diag(D1));
exp_D1_sqrt = diag(exp(pi*1i*d1_sqrt/rho));
phi1 = dft(:, 2);
exp_D1_sqrt_phi1 = exp_D1_sqrt*phi1;
translation = idft*exp_D1_sqrt_phi1;

% The graph heat kernel is the matrix exponential exp(-tL)
g_01 = expm(-.01*L);    
g_1 = expm(-.1*L);
g = expm(-L);

% The GSTFT is defined as exp(-tL)*(V1^*)
gft_01 = g_01 * V1;
gft_1 = g_1 * V1;
gft = g_1 * V1;

gft_01_dft = g_01 * dft;
gft_1_dft = g_1 * dft;
gft_dft = g_1 * dft;
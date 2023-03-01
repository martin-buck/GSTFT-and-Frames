clear
clc

% Circle Graph Laplacian
A = [0 1; 1 0];
D = diag([1 1]);
L = D - A;

% Calculate the eigenvectors
[V1, D1] = eig(L);
v1 = V1(:, 1);
v2 = V1(:, 2);
v1_c = conj(v1)';
v2_c = conj(v2)';

% DFT and IDFT matrices filled with 5th roots of unity and their powers
dft = dftmtx(2);    
idft = dft';

% The graph heat kernel is the matrix exponential exp(-tL) for three values
% of time 't'

g_t = @(t) expm(-t*L);

% Create the analysis operator
eig_diff_arr = [];
count = 1;
for t=0:.1:10
    g = g_t(t);
    A = [g(1,1)*v1_c; g(1,2)*v2_c; g(2,1)*v1_c; g(2,2)*v2_c];
    S = A*ctranspose(A);
    D = eig(S);
    eig_diff = max(D)-min(D);
    eig_diff_arr(count) = eig_diff;
    count = count + 1;
end
plot(eig_diff_arr)
xlabel('t')
ylabel('Maximum Eigenvalue Difference')
title('Eigenvalue Difference of Heat Kernel Frame Operator, N=2')

%% Circle(N)

% Now repeat the above analysis but for circle graphs of varying sizes.
% Create the graph Laplacian for size N

tol = 1e-12;
for i=[5,6,7,8,9,10,50]
    s = [];
    t = [];
    for j=1:i
        if j ~= i
            s(j)=j;
            t(j)=j+1;
        end
    end
    s(j)=j;
    t(j)=1;
    
    G = graph(s, t);
    L = full(laplacian(G));
    [V1, D1] = eig(L);
    V1_conj = ctranspose(V1);
    
    g_t = @(t) expm(-t*L);

    % Create the analysis and frame operator
    eig_diff_arr = [];
    count = 1;
    for t=0:1:20
        g = g_t(t);
        g_T = g';
        A = zeros(i^2, i);
        for j=1:i^2
            if j <= i
                A(j, :) = g_T(j)*V1_conj(j, :);
            elseif j > i && mod(j, i) ~= 0
                A(j, :) = g_T(j)*V1_conj(mod(j, i), :);
            elseif j > i && mod(j, i) == 0
                A(j, :) = g_T(j)*V1_conj(mod(j, i)+1, :);
            end
        end
        S = A*ctranspose(A);
        D = eig(S);
        eig_diff = max(D)-min(D(D>tol));
        eig_diff_arr(count) = eig_diff;
        count = count + 1;
    end
    plot(eig_diff_arr, 'o-')
    xlabel('t')
    ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|')
    legend('N=5', 'N=6', 'N=7', 'N=8', 'N=9', 'N=10', 'N=50')
    title('Frame Operator Eigenvalue Gap')
    grid on
    hold on
end





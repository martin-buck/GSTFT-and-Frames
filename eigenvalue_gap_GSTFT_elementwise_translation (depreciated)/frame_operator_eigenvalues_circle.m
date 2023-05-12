clear
clc
close all

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
for t=0:1:10
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
clear
clc
close all

% Now repeat the above analysis but for circle graphs of varying sizes.
% Create the graph Laplacian for size N

tol = 1e-12;
for i=[5,6,7,8,9,10]
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
    a = 1/sqrt(i)*ctranspose(dftmtx(i));
    
    g_t = @(t) expm(-t*L);

    % Create the analysis and frame operator
    eig_diff_arr = [];
    count = 1;
    times = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
    times_str = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};
    for t = times
        g = g_t(t);
        g_T = g';
        
        % plot the graph heat kernel
        if i == 10
        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(g)
        colorbar;
        title(strcat('Graph Heat Kernel, t=', times_str(count)), 'Fontsize', 36, 'Interpreter', 'latex');
        % saveas(gcf, strcat('graph_heat_kernel_circle', ', t=', times_str{count}, '.png'))
        end

        A = zeros(i^2, i);
        for j=1:i^2
            if j <= i
                A(j, :) = g_T(j)*V1_conj(j, :);
            elseif j > i && mod(j, i) ~= 0
                A(j, :) = g_T(j)*V1_conj(mod(j, i), :);
            elseif j > i && mod(j, i) == 0
                A(j, :) = g_T(j)*V1_conj(mod(j, i)+i, :);
            end
        end
        S = A*ctranspose(A);
        T = ctranspose(A)*A;
        D = eig(S);
        eig_diff = max(D)-min(D(D>tol));
        eig_diff_arr(count) = eig_diff;
        count = count + 1;
    end
    
    plot(eig_diff_arr, 'o-', 'LineWidth', 2)
    xlabel('t','Fontsize', 18)
    xticklabels(times_str)
    ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|', 'Fontsize', 18)
    legend('N=5', 'N=6', 'N=7', 'N=8', 'N=9', 'N=10', 'N=50')
    title('Frame Operator Eigenvalue Gap', 'Interpreter', 'latex', 'Fontsize', 36)
    grid on
    hold on
end

%% Outer-Product Formulation

% Make sure that the above frame matrix A is formed correctly. It SHOULD be
% a multiple of the identity for the circle graph on N vertices

B = zeros(i);

for j=1:i
    B = B + (norm(g(:,j))^2) * V1(:,j) * V1_conj(j,:);
end
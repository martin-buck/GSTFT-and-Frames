3%% Peterson Graph
clear
clc
close all

i = 10;
tol = 1e-12;

A = [0, 0, 1, 1, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 0, 1, 0, 0, 0
    1, 0, 0, 0, 1, 0, 0, 1, 0, 0
    1, 1, 0, 0, 0, 0, 0, 0, 1, 0
    0, 1, 1, 0, 0, 0, 0, 0, 0, 1
    1, 0, 0, 0, 0, 0, 1, 0, 0, 1
    0, 1, 0, 0, 0, 1, 0, 1, 0, 0
    0, 0, 1, 0, 0, 0, 1, 0, 1, 0
    0, 0, 0, 1, 0, 0, 0, 1, 0, 1
    0, 0, 0, 0, 1, 1, 0, 0, 1, 0];

G = graph(A);
plot(G)
L = full(laplacian(G));
[V, D] = eig(L);

V1_conj = ctranspose(V);
g_t = @(t) expm(-t*L);

eig_diff_arr = [];
count = 1;
times = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
times_str = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};

for t = times
    g = g_t(t);
    
    % plot the graph heat kernel
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(g)
    colorbar;
    title(strcat('Graph Heat Kernel, t=', times_str(count)));
    saveas(gcf, strcat('graph_heat_kernel_peterson', ', t=', times_str{count}, '.png'))

    g_T = g';
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

plot(eig_diff_arr, 'o-')
xlabel('t')
xticklabels(times_str)
ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|')
title('Frame Operator Eigenvalue Gap')
grid on
hold on
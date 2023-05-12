% This is now correct and concise code to compute the eigenvalues of the
% GSTFT frame operator/matrix

clear
clc
close all

i = 10;
tol = 1e-12;

plot_kernel = 0;

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
[V, U] = eig(L);

g_t = @(t) expm(-t*L);

count = 1;
times = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
times_str = {'1e-10', '1e-9', '1e-8', '1e-7', '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};
eig_diff_arr = zeros(length(times), 1);

for t = times
    g = g_t(t);
    
    % plot the graph heat kernel
    if plot_kernel == 1
        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(g)
        colorbar;
        title(strcat('Graph Heat Kernel, t=', times_str(count)));
        saveas(gcf, strcat('graph_heat_kernel_peterson', ', t=', times_str{count}, '.png'))
    end

    D = vecnorm(g).^2;
    eig_diff = max(D)-min(D(D>tol));
    eig_diff_arr(count) = eig_diff;
    count = count + 1;
end

plot(eig_diff_arr, 'o-', 'LineWidth', 2)
xlabel('t','Fontsize', 18)
xticklabels(times_str)
ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|', 'Fontsize', 18)
legend('N=10', 'N=100', 'N=1000')
title('Frame Operator Eigenvalue Gap', 'Interpreter', 'latex', 'Fontsize', 36)
grid on
% This is now correct and concise code to compute the eigenvalues of the
% GSTFT frame operator/matrix

clear
clc
close all

plot_kernel = 1;

% load Minnisota road graph
M = load('road-minnesota/minnesota.mat');

% now create a graph from the node-node pairs in the above list we imported
A = M.Problem.A;
coord = M.Problem.aux.coord;
G = graph(A);
G = rmnode(G, [348, 349]);
% these two nodes are their own connected component, apparently?
coord(348, :) = [];
coord(348, :) = [];
figure('units','normalized','outerposition',[0 0 1 1])
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'NodeLabel', 1:2640)

tol = 1e-12;

A = full(adjacency(G));
D = diag(sum(A, 2));
L = D-A;
L = full(L);
N = length(A(:,1));

g_t = @(t) expm(-t*L);

count = 1;
times = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
times_str = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};
eig_diff_arr = zeros(length(times), 1);

for t = times
    g = g_t(t);
    g_T = g';
    
    if plot_kernel == 1
        % plot the graph heat kernel
        figure('units', 'normalized', 'outerposition', [0 0 1 1])
        imagesc(g)
        colorbar;
        title(strcat('Graph Heat Kernel, t=', times_str(count)));
        saveas(gcf, strcat('graph_heat_kernel_minn', ', t=', times_str{count}, '.png'))
    end

    D = vecnorm(g_T).^2;
    eig_diff = max(D)-min(D(D>tol));
    eig_diff_arr(count) = eig_diff;
    count = count + 1;
end

figure
plot(eig_diff_arr, 'o-', 'LineWidth', 2)
xlabel('t','Fontsize', 18)
xticklabels(times_str)
ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|', 'Fontsize', 18)
title('Frame Operator Eigenvalue Gap', 'Interpreter', 'latex', 'Fontsize', 36)
grid on
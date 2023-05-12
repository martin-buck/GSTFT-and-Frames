clear
clc
close all

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

[V1, W1] = eig(full(L));

V1_conj = ctranspose(V1);
g_t = @(t) expm(-t*L);

% At this point we need to create the analysis or frame operator which
% requires creating an N^2 by N matrix which is too big to store locally!
% How to we get around this? Use sparse matrices and write as a sum of
% outer products so we don't have to form the full N^2 by N matrix

eig_diff_arr = [];
count = 1;
times = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
times_str = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};
for t = times
    g = g_t(t);
    g_T = g';
    
    % plot the graph heat kernel
    figure('units', 'normalized', 'outerposition', [0 0 1 1])
    imagesc(g)
    colorbar;
    title(strcat('Graph Heat Kernel, t=', times_str(count)));
    saveas(gcf, strcat('graph_heat_kernel_path', ', t=', times_str{count}, '.png'))
    
    % outer product formulation
    T = zeros(N, N);
    for j=1:N^2
        if j <= N
            rowT = g_T(j)*V1_conj(j, :);
        elseif j > N && mod(j, N) ~= 0
            rowT = g_T(j)*V1_conj(mod(j, N), :);
        elseif j > N && mod(j, N) == 0
            rowT = g_T(j)*V1_conj(mod(j, N)+N, :);
        end
        colT_conj = ctranspose(rowT);
        T = T + colT_conj*rowT;
        
        if mod(j, 100) == 0
            disp(j)
        end
    end
    D = eig(T);
    eig_diff = max(D)-min(D(D>tol));
    eig_diff_arr(count) = eig_diff;
    count = count + 1;
end
close all

plot(eig_diff_arr, 'o-')
xlabel('t')
xticklabels(times_str)
ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|')
title('Frame Operator Eigenvalue Gap')
saveas(gcf, strcat('maximum_eig_diff', ', t=', times_str{count}, '.png'))
grid on
hold on
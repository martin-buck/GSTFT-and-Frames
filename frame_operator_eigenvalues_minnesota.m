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
coord(348, :) = [];
coord(348, :) = [];
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
% How to we get around this?
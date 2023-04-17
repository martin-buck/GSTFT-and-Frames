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

%% Shuman/Vanderghenyst translation

% Now build the graph Laplacian and get its eigenvectors to create the
% fourier exponential function that is used in Shuman et al (2013) to
% demonstrate translation in their paper on this graph
A = full(adjacency(G));
D = diag(sum(A, 2));
L = D-A;
L = full(L);
N = length(A(:,1));

[V, W] = eig(full(L));

% V* is our FT matrix
FT = ctranspose(V);
IFT = V;

%f_FT = exp(-5*diag(W));
%f_FT_1 = f_FT/norm(f_FT);
%f = IFT*f_FT_1;
f = zeros(N,1); f(1)=1;
f_FT = FT*f;
f_FT_1 = f_FT/norm(f_FT);

% plot the function
figure
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', f, 'NodeCData', f)
colormap default
colorbar
title('f')

% plot the FT of the function
figure
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', f_FT, 'NodeCData', f_FT)
colormap default
colorbar
title('FT(f)')

figure
plot(f_FT, 'o-')
grid on
xlabel('i')
ylabel('\lambda_i')
title('Laplacian Spectrum')

% they translate T_200, T_1000, and T_2000
f_t200_conv = f_FT_1.*FT(:, 200);
f_t200 = sqrt(N)*IFT*f_t200_conv;
figure
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', f_t200, 'NodeCData', f_t200)
colormap default
colorbar
title('T_{200}f')

f_t1000_conv = f_FT_1.*FT(:, 1000);
f_t1000 = sqrt(N)*IFT*f_t1000_conv;
figure
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', f_t1000, 'NodeCData', f_t1000)
colormap default
colorbar
title('T_{1000}f')

f_t2000_conv = f_FT_1.*FT(:, 2000);
f_t2000 = sqrt(N)*IFT*f_t2000_conv;
figure
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', f_t2000, 'NodeCData', f_t2000)
colormap default
colorbar
title('T_{2000}f')

f_t10_conv = f_FT_1.*FT(:, 10);
f_t10 = sqrt(N)*IFT*f_t10_conv;
figure
plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', f_t10, 'NodeCData', f_t10)
colormap default
colorbar
title('T_{10}f')

close all
%% GSTFT/Frame Translation in Fourier Domain (NxN block)

% graph heat kernel. Cannot naively create an N^2 by N matrix for this
% graph as this exceeds memory requirements
N_max = 10;
g_t = @(t) expm(-t*full(L));

% now create an N^2 by N matrix that consists of N by N blocks that are
% scaled versions of the Fourier transform matrix with scalings coming from
% the columns of the graph heat kernel g_t
t = .1;
g = g_t(t);

% the way f is defined in Shuman et al. (2013) is via its fourier transform
for i=201:201+N_max+1
    gstft_block = diag(g(:, i-1))*FT;
    % this isn't the inverse gstft as defined in Overleaf but just the
    % inverse fourier transform of this particular block
    igstft_f = sqrt(N)*IFT*gstft_block*f;
    
    % what does the gstft look like on this particular block without
    % inverse transforming back?
    gstft_f = gstft_block*f;
    
    % plot igstft 3D
    figure
    plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', igstft_f, 'NodeCData', igstft_f)
    colormap default
    colorbar
    set(gcf, 'Position', get(0, 'Screensize'));
    title_str = strcat('GSTFT ', ' Block ', num2str(i-1), ', t=', 'one-tenth');
    title(title_str);
    saveas(gcf, title_str, 'png');
    
    % plot igstft flat
    figure
    plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'NodeCData', igstft_f)
    colormap default
    colorbar
    set(gcf, 'Position', get(0, 'Screensize'));
    title_str = strcat('GSTFT ', ' Block ', num2str(i-1), ', t=', 'one-tenth', ', Flat');
    title(title_str);
    saveas(gcf, title_str, 'png');
    
    % plot gstft 3D
    figure
    plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'Zdata', gstft_f, 'NodeCData', gstft_f)
    colormap default
    colorbar
    set(gcf, 'Position', get(0, 'Screensize'));
    title_str = strcat('GSTFT Fourier Domain', ' Block ', num2str(i-1), ', t=', 'one-tenth');
    title(title_str);
    saveas(gcf, title_str, 'png');
    
    % plot gstft flat
    figure
    plot(G, 'Xdata', coord(:,1), 'Ydata', coord(:,2), 'NodeCData', gstft_f)
    colormap default
    colorbar
    set(gcf, 'Position', get(0, 'Screensize'));
    title_str = strcat('GSTFT Fourier Domain', ' Block ', num2str(i-1), ', t=', 'one-tenth');
    title(title_str);
    saveas(gcf, title_str, 'png');
end
close all
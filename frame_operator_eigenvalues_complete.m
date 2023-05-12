% This is now correct and concise code to compute the eigenvalues of the
% GSTFT frame operator/matrix

clear
clc
close all

plot_kernel = 0;

% Now repeat the above analysis but for circle graphs of varying sizes.
% Create the graph Laplacian for size N
tol = 1e-12;
for i=[10, 100, 1000]
    v = zeros(i, 1);
    for j=1:i
        v(j) = i;
    end
    L = -1*ones(i);
    L = L + diag(v);
    [V, U] = eig(L);
    
    g_t = @(t) expm(-t*L);

    % Create the analysis and frame operator
    count = 1;
    times = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
    times_str = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};
    eig_diff_arr = zeros(length(times), 1);
    for t = times
        g = g_t(t);
        g_T = g';
        
        if plot_kernel == 1
            % plot the graph heat kernel
            figure('units','normalized','outerposition',[0 0 1 1])
            imagesc(g)
            colorbar;
            title(strcat('Graph Heat Kernel, t=', times_str(count)));
            saveas(gcf, strcat('graph_heat_kernel_complete', ', t=', times_str{count}, '.png'))
        end
        
        D = vecnorm(g_T).^2;
        eig_diff = max(D)-min(D(D>tol));
        eig_diff_arr(count) = eig_diff;
        count = count + 1;
    end
    
    plot(eig_diff_arr, 'o-', 'LineWidth', 2)
    xlabel('t', 'Fontsize', 18)
    xticklabels(times_str)
    ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|', 'Fontsize', 18)
    legend('N=10', 'N=100', 'N=1000')
    title('Frame Operator Eigenvalue Gap', 'FontSize', 36, 'Interpreter', 'latex')
    grid on
    hold on
end
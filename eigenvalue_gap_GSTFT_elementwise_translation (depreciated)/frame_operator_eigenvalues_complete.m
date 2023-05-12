%% Complete(N)
clear
clc
close all

% Now repeat the above analysis but for circle graphs of varying sizes.
% Create the graph Laplacian for size N
tol = 1e-12;
for i=[5, 6, 7, 8, 9, 10]
    v = zeros(i, 1);
    for j=1:i
        v(j) = i;
    end
    L = -1*ones(i);
    L = L + diag(v);
    [V1, D1] = eig(L);
    
    V1_conj = ctranspose(V1);
    a = 1/sqrt(i)*  ctranspose(dftmtx(i));
    
    g_t = @(t) expm(-t*L);

    % Create the analysis and frame operator
    eig_diff_arr = [];
    count = 1;
    times = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000];
    times_str = {'1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1', '10', '100', '1000', '10000', '100000'};
    for t = times
        g = g_t(t);
        g_T = g';
        
%         % plot the graph heat kernel
%         figure('units','normalized','outerposition',[0 0 1 1])
%         imagesc(g)
%         colorbar;
%         title(strcat('Graph Heat Kernel, t=', times_str(count)));
%         saveas(gcf, strcat('graph_heat_kernel_complete', ', t=', times_str{count}, '.png'))

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
    xlabel('t', 'Fontsize', 18)
    xticklabels(times_str)
    ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|', 'Fontsize', 18)
    legend('N=5', 'N=6', 'N=7', 'N=8', 'N=9', 'N=10', 'N=50')
    title('Frame Operator Eigenvalue Gap', 'FontSize', 36, 'Interpreter', 'latex')
    grid on
    hold on
end
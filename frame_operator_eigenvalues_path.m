%% Path(N)
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
    
    G = graph(s, t);
    L = full(laplacian(G));
    [V1, D1] = eig(L);
    
    V1_conj = ctranspose(V1);
    a = 1/sqrt(i)*  ctranspose(dftmtx(i));
    
    g_t = @(t) expm(-t*L);

    % Create the analysis and frame operator
    eig_diff_arr = [];
    count = 1;
    t_arr = [0:.1:10];
    for t=t_arr
        g = g_t(t);
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
    plot(t_arr, eig_diff_arr, 'o-')
    xlabel('t')
    ylabel('Maximum Eigenvalue Difference |\lambda_{max}-\lambda_{min}|')
    legend('N=5', 'N=6', 'N=7', 'N=8', 'N=9', 'N=10', 'N=50')
    title('Frame Operator Eigenvalue Gap')
    grid on
    hold on
end
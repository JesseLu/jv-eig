%% Rayleigh Quotient Iteration 
% An algorithm from page 207 of Numerical Linear Algebra, Trefethen and Bau.

function [v, l] = rqi(A, sAinv, v, max_iters, err_lim)
    v = v / norm(v);
    l = v' * A(v);

    for k = 1 : max_iters+1
        err(k) = norm(A(v) - l * v); % Compute error.

        semilogy(0:(k-1), err, '.-'); % Plot error.
        ylabel('Eigenvector error'); xlabel('iteration'); drawnow; 

        if (err(k) < err_lim) || (k >= max_iters) % Check if we're done
            break
        end

        w = sAinv(l, v); % Solve for new eigenvector guess (inverse iteration).
        v = w / norm(w); % Normalize.
        l = v' * A(v); % Solve for new eigenvalue guess (rayleigh quotient).
    end
    
        

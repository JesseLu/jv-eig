
function [x, err] = qmr(A, b, max_iters)

    m_dot = @(x, y) dot(conj(x), y);
    n = length(b);
    x = zeros(n, 1);
    x_Q = x;
 %    xhat = zeros(n, 1);
    r = b - A * x;
    rhat = r;
    p = r;
    phat = rhat;
    rho = zeros(max_iters, 1);
    rho(1) = m_dot(rhat, r);
    err = zeros(max_iters, 1);

    % Extra variables needed for QMR.
    q = zeros(n,1);
    tau = norm(r);
    theta = zeros(max_iters, 1);

    term_err = 1e-6 * norm(b);

    for k = 1 : max_iters

        err(k) = norm(r);
        fprintf('%d: %e\n', k, err(k));
        if mod(k, 400) == 0
            subplot 211; plot([real(x), imag(x), abs(x)], '.-');
            subplot 212; plot([real(x_Q), imag(x_Q), abs(x_Q)], '.-'); drawnow;
        end

        if err(k) < term_err
            err = err(1:k);
            return
        end

        % rho(k) = m_dot(rhat, r);

        % Step 1.
        v = A * p;
        vhat = A.' * phat;

        alpha = rho(k) / m_dot(phat, v);
        x = x + alpha * p;
%         xhat = xhat + alpha * phat;

        r = r - alpha * v;
        rhat = rhat - alpha * vhat;


        % Step 2 (QMR).
        theta(k+1) = norm(r) / tau;
        c = 1 / sqrt(1 + theta(k+1)^2);
        tau = tau * theta(k+1) * c;
        q = c^2 * theta(k)^2 * q + c^2 * alpha * p;
        x_Q = x_Q + q;



        % Step 3.
        rho(k+1) = m_dot(rhat, r);
        beta = rho(k+1) / rho(k);

        p = r + beta * p;
        phat = rhat + beta * phat;
    end
    err = err(1:k);


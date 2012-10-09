function [omega, E, H, err] = eigenmode(sim, omega, E, ...
                                    d_prim, d_dual, s_prim, s_dual, ...
                                    mu, epsilon,  ...
                                    max_iters, err_lim, varargin);

    %% Transform problem into eigenvalue problem.
    % The eigenvalue matrix is given by:
    %     1/sqrt(epsilon) curl curl 1/sqrt(epsilon)
    % and this is from the transformation F = sqrt(epsilon) E.

    % Matrices needed for the transform.
    [A1, A2, m, e, b] = fds_matrices(omega, s_prim, s_dual, mu, epsilon, mu);

    dims = size(E{1});
    my_zeros = {zeros(dims), zeros(dims), zeros(dims)};
    
    % Multiplication be the eigenvalue matrix.
    multA = @(v) (e.^-0.5) .* (A1 * (A2 * (e.^-0.5 .* v)));

    % Computes (A - lambda I)^-1 v using Maxwell.
    sAinv = @(l, v) shifted_A_inverse(l, v, e, d_prim, d_dual, s_prim, s_dual, ...
                                        mu, epsilon, my_zeros, sim);
    sAinv_err = @(l, v, w) norm(multA(w) - l * w - v); % Measure of SAinv error.

    lambda = omega^2;
    v = sqrt(e) .* [E{1}(:); E{2}(:); E{3}(:)];
    
    if isempty(varargin)
        [v, lambda] = rqi(multA, sAinv, v, max_iters, err_lim);
    elseif strcmp(varargin{1}, 'inviter')
        [v, lambda] = inviter(multA, sAinv, v, lambda, max_iters, err_lim);
    elsle
        error('Invalid optional parameter.');
    end


    % Get the frequency of the mode.
    omega = sqrt(lambda);

    % Transform back to E-field.
    x = v ./ sqrt(e);
    N = prod(dims);
    E = {reshape(x(1:N), dims), ...
        reshape(x(N+1:2*N), dims), ...
        reshape(x(2*N+1:3*N), dims)};

    % Calculate H-field from E-field.
    y = A2 * x ./ (-i * omega * m);
    H = {reshape(y(1:N), dims), ...
        reshape(y(N+1:2*N), dims), ...
        reshape(y(2*N+1:3*N), dims)};

    % Calculate errors.
    err.actual = sAinv_err(lambda, 0, v); % Final error.
    err.E = norm(A1 * ((1./m) .* (A2 * x)) - omega^2 * e .* x) / norm(x); 
    err.H = norm(A2 * ((1./e) .* (A1 * y)) - omega^2 * m .* y) / norm(y);

function [x] = shifted_A_inverse(shift, v, e, d_prim, d_dual, s_prim, s_dual, ...
                                    mu, epsilon, E0, sim)
    dims = size(E0{1});
    N = prod(dims);
    b = (sqrt(e) .* v) ./ (-i * sqrt(shift));
    J = {reshape(b(1:N), dims), ...
            reshape(b(N+1:2*N), dims), ...
            reshape(b(2*N+1:3*N), dims)};
%     % Override.
%     v = v ./ sqrt(e);
%     E0 = {reshape(v(1:N), dims), ...
%             reshape(v(N+1:2*N), dims), ...
%             reshape(v(2*N+1:3*N), dims)};
    [E, H, err] = sim(sqrt(shift), J);
    x = sqrt(e) .* [E{1}(:); E{2}(:); E{3}(:)];
    

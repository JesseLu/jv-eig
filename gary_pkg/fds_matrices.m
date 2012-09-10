%% fds_matrices
% Advanced feature which constructs the matrices and vectors that make up the simulation. 

%% Syntax
% [A1, A2, mu, epsilon, b] = fds_matrices(omega, s_prim, s_dual, mu, epsilon, J)

%% Description
% Produces the matrices and vectors needed to form the numerical model locally.
% This is an advanced feature that is useful if we want to locally check the
% error of the solution, or compute additional problem information such as
% derivative information.

%% Input parameters.
% Input parameters exactly match those used in the |fds.m| function.

%% Output parameters.
% * |A1|, |A2| -- matrices (sparse).
% * |m|, |e|, |b| -- vectors.

% Together, these allow one to form the problem on one's own machine.

%% Example: Local error check.
% Using the same parameters as from an |fds()| call, we can then compute the
% error in the simulation locally. Note that the locally computed error may
% differ from the server-computed value for implementation reasons.

%   % Run the simulation.
%   [E, H, err] = fds(omega, s_prim, s_dual, mu, epsilon, E, J, max_iters, err_thresh, 'text');
%
%   x = ([E{1}(:) ; E{2}(:) ; E{3}(:)]); % "Vectorize" E-field values.
%   y = ([H{1}(:) ; H{2}(:) ; H{3}(:)]); % "Vectorize" H-field values.
%
%   % Construct matrices
%   [A1, A2, m, e, b] = fds_matrices(omega, s_prim, s_dual, mu, epsilon, J); 
%
%   % Calculate errors in the wave equation for both E-fields and H-fields.
%   norm(A1 * ((1./m) .* (A2 * x)) - omega^2 * e .* x - b) / norm(b) % Error for E-field wave equation.
%   norm(A2 * ((1./e) .* (A1 * y)) - omega^2 * m .* y - A2 * (b ./ (-i*omega*e))) / norm(b) % Error for H-field wave equation.

function [A1, A2, m, e, b] = fds_matrices(omega, s_prim, s_dual, mu, epsilon, J)

    %% Get the shape.
    dims = size(epsilon{1});
    N = prod(dims);
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    my_blkdiag = @(z) blkdiag(my_diag(z{1}), my_diag(z{2}), my_diag(z{3}));

    %% Get the relevant derivative matrices.
    [spx, spy, spz] = ndgrid(s_prim{1}, s_prim{2}, s_prim{3});
    [sdx, sdy, sdz] = ndgrid(s_dual{1}, s_dual{2}, s_dual{3});

    Dx = deriv('x', dims); % Derivative in x, y, and z directions.
    Dy = deriv('y', dims);
    Dz = deriv('z', dims);
    Z = sparse(N, N);

    % Forward differences (used to compute H from E).
    Dfx = my_diag(sdx.^-1) * Dx;
    Dfy = my_diag(sdy.^-1) * Dy;
    Dfz = my_diag(sdz.^-1) * Dz;

    % Backward differences (used to compute E from H).
    Dbx = -my_diag(spx.^-1) * Dx';
    Dby = -my_diag(spy.^-1) * Dy';
    Dbz = -my_diag(spz.^-1) * Dz';

    %% Form matrices
    A1 = [  Z, -Dbz, Dby; ...
            Dbz, Z, -Dbx; ...
            -Dby, Dbx, Z];

    A2 = [  Z, -Dfz, Dfy; ...
            Dfz, Z, -Dfx; ...
            -Dfy, Dfx, Z];

    m = [mu{1}(:) ; mu{2}(:) ; mu{3}(:)];
    e = [epsilon{1}(:) ; epsilon{2}(:) ; epsilon{3}(:)];


    b = -i * omega * [J{1}(:) ; J{2}(:) ; J{3}(:)];



%% Derivative for three dimensional space.
% Note that we are making the forward derivative only.
% Also, we assume periodic boundary conditions.
function [D] = deriv(dir, shape)

shift = (dir == 'xyz'); % Direction of shift.

% Get the displaced spatial markers.
my_disp = @(n, shift) mod([1:n] + shift - 1, n) + 1;
[i, j, k] = ndgrid(my_disp(shape(1), shift(1)), ...
                    my_disp(shape(2), shift(2)), ...
                    my_disp(shape(3), shift(3)));

% Translate spatial indices into matrix indices.
N = prod(shape);
i_ind = 1 : N;
j_ind = i + (j-1) * shape(1) + (k-1) * shape(1) * shape(2);

% Create the sparse matrix.
D = sparse([i_ind(:); i_ind(:)], [i_ind(:), j_ind(:)], ...
            [-ones(N,1); ones(N,1)], N, N);

function [omega, E, H, err] = arm_eig(type, E0, amp, ds, base_omega, max_iters, err_lim)
    
    path(path, genpath('../fds-client/'));

    v = load('armand1.mat');
    epsilon = v.epsilon;

    % Get the initial E-field (and other values) needed to find the eigenmode.
    [omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E, J, sim] = ...
        get_case(type, E0, amp, ds, base_omega, double(epsilon));
        
    % Find the eigenmode.
    [omega, E, H, err] = eigenmode(sim, omega, E, ...
                                    d_prim, d_dual, s_prim, s_dual, ...
                                    mu, epsilon, ...
                                    max_iters, err_lim);


    fprintf('Omega: %1.3e + i%1.3e\nErrors: %e (actual), %e (E), %e (H)\n', ...
            real(omega), imag(omega), err.actual, err.E, err.H);

function [omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E, J, sim_eig] = ...
            get_case(type, E0, amp, ds, base_omega, epsilon, varargin)

    %% Form the excitation source.
    spr = [-1:1];
    pos = {[352; 158], [455; 194], [350; 230], [245; 196]};
    if strcmp(type, 'iso')
        cnt = find(amp == 1);
        cnt = cnt(1); 
        epsilon = epsilon(pos{cnt}(1)+[-75:75], pos{cnt}(2)+[-36:36], :);
        dims = size(epsilon);
        J = {zeros(dims), zeros(dims), zeros(dims)}; 
        c = round(dims/2);
        J{2}(c(1)+spr, c(2)+spr, c(3)+spr) = 1;
    elseif strcmp(type, 'point')
        dims = size(epsilon)
        J = {zeros(dims), zeros(dims), zeros(dims)}; 
        for k = 1 : length(pos) 
                % pos{k} = pos{k} - 100;
                J{2}(pos{k}(1)+spr, pos{k}(2)+spr, dims(3)/2+spr) = amp(k);
        end
    elseif strcmp(type, 'initialE')
        dims = size(epsilon);
        J = {zeros(dims), zeros(dims), zeros(dims)}; 
    else 
        error('Invalid TYPE.');
    end

    %% Downsample both the frequency, epsilon, and excitation.
    omega = base_omega * ds;
    for k = 1 : 3
        J{k} = J{k}(ds:ds:end, ds:ds:end, ds:ds:end);
    end
    epsilon = epsilon(ds:ds:end, ds:ds:end, ds:ds:end);
    epsilon = eps_smooth(epsilon);
    dims = size(epsilon{1});
    
    %% Form the other necessary operators for fdfd simulation.
    d = {ones(dims(1), 1), ones(dims(2), 1), ones(dims(3), 1)};
    d_prim = d;
    d_dual = d;
    [s_prim, s_dual] = make_scpml(omega, dims, 10/ds);
    mu = {ones(dims), ones(dims), ones(dims)}; 
    E = {zeros(dims), zeros(dims), zeros(dims)}; 

    %% Form the function handle to run the simulation.
    sim_fdfd = @(omega, J, E) my_sim(omega, d_prim, d_dual, s_prim, s_dual, ...
                        mu, epsilon, E, J, ...
                        1e5, 1e-6, 'plot', ...
                        2, round(dims(3)/2), ['amp ', num2str(amp(:)'), ' ds ', num2str(ds)]);

    %% The function handle that will be passed to the eigenmode function.
    sim_eig = @(omega, J) sim_fdfd(omega, J, E);

    %% Debugging.
    subplot 321;
    zc = round(dims(3)/2);
    imagesc(epsilon{3}(:,:,zc)'); axis equal tight;
    imagesc(((1 - J{2}(:,:,zc)) .* epsilon{2}(:,:,zc))'); axis equal tight;
    % pause

    %% Get guess field.
    if strcmp(type, 'initialE')
        for k = 1 : 4
            c = round(pos{k}/ds); % Center the things here.
            c(3) = round(dims(3)/2);
            for l = 1 : 3
                s = size(E0{k}{l});
                bi = round(c(:)-s(:)/2); % start indices.
                ei = round(c(:)+s(:)/2)-1; % end indices.
                E{l}(bi(1):ei(1), bi(2):ei(2), bi(3):ei(3)) = ...
                    E{l}(bi(1):ei(1), bi(2):ei(2), bi(3):ei(3)) + amp(k) * E0{k}{l};
            end
        end
            
    else
        title('Initial simulation');
        [E, H, err] = sim_eig(omega, J);
%         if strcmp(type, 'iso')
%             for k = 1 : 3
%                 D{k} = zeros(dims);
%                 D{k}(
    end 
    % subplot 111
    imagesc(E{2}(:,:,zc)'); axis equal tight;
    % pause
    
%% Calculates s_prim and s_dual for a regularly spaced grid of dimension DIMS.
% Grid spacing assumed to be 1.
function [s_prim, s_dual] = make_scpml(omega, dims, t_pml)
    % Helper functions.
    pos = @(z) (z > 0) .* z; % Only take positive values.
    l = @(u, n) pos(t_pml - u) + pos(u - (n - t_pml)); % Distance to nearest pml boundary.

    % Compute the stretched-coordinate grid spacing values.
    for k = 1 : 3
        s_prim{k} = 1 - i * (4 / omega) * (l(0:dims(k)-1, dims(k)) / t_pml).^4;
        s_dual{k} = 1 - i * (4 / omega) * (l(0.5:dims(k)-0.5, dims(k)) / t_pml).^4;
    end


function [E, H, err] = my_sim(omega, d_prim, d_dual, s_prim, s_dual, ...
                                mu, epsilon, E0, J, ...
                                fds_max_iters, fds_err_lim, view_option, ...
                                E_index, E_loc, name)
    subplot 322;
    [E, H, err] = fds(omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E0, J, ...
                        fds_max_iters, fds_err_lim, view_option);

    subplot(3,1,2:3);
    imagesc(abs(E{E_index}(:,:,E_loc)')); axis equal tight;

    subplot 321;
    drawnow
    saveas(gcf, [name, ' ', datestr(now, 'mm-dd-HH:MM:SS')], 'png')

function [epsilon] = eps_smooth(e)
    epsilon{1} = 0.5 * (e + e([2:end, 1], :, :));
    epsilon{2} = 0.5 * (e + e(:, [2:end, 1], :));
    epsilon{3} = 0.5 * (e + e(:, :, [2:end, 1]));
    

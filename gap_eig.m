function [omega, E, H, err] = gap_eig(case_name, exc_name, max_iters, err_lim, varargin)
    path(path, genpath('../fds-client/'));
    if ~isempty(varargin)
        [omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E, J, sim] = ...
            get_case(case_name, exc_name, varargin{1});
    else
        [omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E, J, sim] = ...
            get_case(case_name, exc_name);
    end
    dims = size(epsilon{1});
    [omega, E, H, err] = eigenmode(sim, omega, E, ...
                                    d_prim, d_dual, s_prim, s_dual, ...
                                    mu, epsilon, ...
                                    max_iters, err_lim);


    fprintf('Omega: %1.3e + i%1.3e\nErrors: %e (actual), %e (E), %e (H)\n', ...
            real(omega), imag(omega), err.actual, err.E, err.H);

function [omega, d_prim, d_dual, s_prim, s_dual, mu, epsilon, E, J, sim2] = ...
            get_case(name, exc_name, varargin)
    base_omega = 0.11;
    switch(name)
        case {'1lo', '1med', '1hi'}
            xyz = 'xyz';
            % Determine downsampling ratio.
            if strcmp(name, '1lo')
                ds = 3; 
            elseif strcmp(name, '1med')
                ds = 2; 
            elseif strcmp(name, '1hi')
                ds = 1; 
            end
            omega = base_omega * ds;
            v = load('sonia1.mat');
            epsilon = {v.ex, v.ey, v.ez};
%             J = v.E;
            for k = 1 : 3
%                 epsilon{k} = permute(hdf5read('L3-epsilon.h5', ['/fields/e', xyz(k)]), ...
%                                         [3 2 1]);
                epsilon{k} = epsilon{k}(ds:ds:end, ds:ds:end, ds:ds:end);
                epsilon{k} = double(epsilon{k});
%                 J{k} = J{k}(ds:ds:end, ds:ds:end, ds:ds:end);
%                 J{k} = double(J{k});
            end
            dims = size(epsilon{1});
            

            d = {ones(dims(1), 1), ones(dims(2), 1), ones(dims(3), 1)};
            d_prim = d;
            d_dual = d;
            [s_prim, s_dual] = make_scpml(omega, dims, 10/ds);
            mu = {ones(dims), ones(dims), ones(dims)}; 

            J = {zeros(dims), zeros(dims), zeros(dims)}; 
            c = ceil(dims/2);
            spr = unique(round([-4/ds:4/ds]));
            J{2}(c(1)+spr, c(2)+spr, c(3)+spr) = 1;

            E = {zeros(dims), zeros(dims), zeros(dims)}; 
            sim = @(omega, J, E) my_sim(omega, d_prim, d_dual, s_prim, s_dual, ...
                                mu, epsilon, E, J, ...
                                1e6, 1e-6, 'plot', ...
                                2, round(dims(3)/2), name);

        case {'2lo', '2med', '2hi'}
            xyz = 'xyz';
            % Determine downsampling ratio.
            if strcmp(name, '2lo')
                ds = 3; 
            elseif strcmp(name, '2med')
                ds = 2; 
            elseif strcmp(name, '2hi')
                ds = 1; 
            end
            omega = 2 * base_omega * ds;
            v = load('sonia2.mat');
            epsilon = {v.ex, v.ey, v.ez};
            J = {0*v.ex, 0*v.ey, 0*v.ez};
            dims = size(J{1});
            cen = round(dims/2);
            for k = -1 : 1
                if strcmp(exc_name, 'Pz0')
                    J{3}(:,:,cen(3)+k) = v.Pz0';
                elseif strcmp(exc_name, 'Pz45')
                    J{3}(:,:,cen(3)+k) = v.Pz45';
                elseif strcmp(exc_name, 'center')
                    J{3}(cen(1)+[-4:4], cen(2)+[-4:4], cen(3)+k) = 1;
                else
                    error('Invalid EXC_NAME.');
                end
            end

            for k = 1 : 3
%                 epsilon{k} = permute(hdf5read('L3-epsilon.h5', ['/fields/e', xyz(k)]), ...
%                                         [3 2 1]);
                epsilon{k} = epsilon{k}(ds:ds:end, ds:ds:end, ds:ds:end);
                epsilon{k} = double(epsilon{k});
                J{k} = J{k}(ds:ds:end, ds:ds:end, ds:ds:end);
%                 J{k} = double(J{k});
            end
            dims = size(epsilon{1});
            

            d = {ones(dims(1), 1), ones(dims(2), 1), ones(dims(3), 1)};
            d_prim = d;
            d_dual = d;
            [s_prim, s_dual] = make_scpml(omega, dims, 10/ds);
            mu = {ones(dims), ones(dims), ones(dims)}; 

%             J = {zeros(dims), zeros(dims), zeros(dims)}; 
%             c = ceil(dims/2);
%             spr = unique(round([-4/ds:4/ds]));
%             J{3}(c(1)+spr, c(2)+spr, c(3)+spr) = 1;

            E = {zeros(dims), zeros(dims), zeros(dims)}; 
            sim = @(omega, J, E) my_sim(omega, d_prim, d_dual, s_prim, s_dual, ...
                                mu, epsilon, E, J, ...
                                1e6, 1e-6, 'plot', ...
                                3, round(dims(3)/2), name);

        otherwise
            error('No case by that name found.');
    end

    sim2 = @(omega, J) sim(omega, J, E);
    subplot 321;
    imagesc(epsilon{3}(:,:,round(dims(3)/2))'); axis equal tight;
%     subplot 111;
    imagesc(J{3}(:,:,round(dims(3)/2))'); axis equal tight;
%     pause
    if ~isempty(varargin)
        E = varargin{1};
    else
        title('Initial simulation');
        [E, H, err] = sim2(omega, J);
    end
%     imagesc(epsilon{1}(:,:,round(dims(3)/2))'); axis equal tight;
    
    
    
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

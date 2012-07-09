function [r] = my_batch(w, params, ds, varargin)
    v = load('armand1.mat');
    epsilon = v.epsilon;

    if isempty(varargin)
        res = arm_indiv(0.076, ds, epsilon)
        save(['my_batch_indiv_', num2str(ds), '.mat'], 'res');
    else
        res = varargin{1};
    end


    E0 = arm_E0(res, round([90 24 30]/ds))
%     imagesc(E0{1}{2}(:,:,25)')
%     pause

%     params = {[1 1 1 1], [1 0 -1 0], [0 1 0 -1], [1 -1 1 -1]};
%     w = [0.076 0.076 0.076 0.0766];
    omega = {};
    E = {};
    H = {};
    err = {};
    eps = {};
    for k = 1 : length(params)
        [omega{k}, E{k}, H{k}, err{k}, eps{k}] = ...
            arm_eig('initialJ', epsilon, E0, params{k}, ds, w(k), 10, 1e-8, [150 100]); 
        save(['my_batch_', num2str(ds), '.mat'], 'params', 'omega', 'E', 'H', 'err', 'eps');
    end
    r.indiv = res;
    r.params = params;
    r.omega = omega;
    r.E = E;
    r.H = H;
    r.eps = eps;
    save(['my_batch_r_', num2str(ds), '.mat'], 'r');

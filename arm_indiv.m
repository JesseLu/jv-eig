% Gets the individual resonances.
function [res] = arm_indiv(w, ds, epsilon)

% v = load('armand1.mat');
% epsilon = v.epsilon;
% ds = 3; % Downsample by a factor of 3.
for k = 1 : 4
    amp = zeros(4, 1);
    amp(k) = 1;
    [omega, E, H, err] = arm_eig('iso', epsilon, [], amp, ds, w(k), 10, 1e-6, [0 0]);
    res{k} = {omega, E, H, err};
end
save quad_script_resd.mat res

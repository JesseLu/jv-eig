% Gets the individual resonances.

v = load('armand1.mat');
epsilon = v.epsilon;
ds = 3; % Downsample by a factor of 3.
for k = 1 : 4
    amp = zeros(4, 1);
    amp(k) = 1;
    [omega, E, H, err] = arm_eig('iso', amp, ds, 0.0766, epsilon, 10, 1e-6);
    res{k} = {omega, E, H, err};
end
save quad_script_resd.mat res

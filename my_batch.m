ds = 3;
v = load('armand1.mat');
epsilon = v.epsilon;

res = arm_indiv(0.076, ds, epsilon)
save batch.mat res

E0 = arm_E0(res, round([90 24 30]/ds))

params = {[1 1 1 1], [1 0 -1 0], [0 1 0 -1], [1 -1 1 -1]};
w = [0.076 0.076 0.076 0.0766];
omega = {};
E = {};
H = {};
err = {};
eps = {};
for k = 1 : length(params)
    [omega{k}, E{k}, H{k}, err{k}, eps{k}] = ...
        arm_eig('initialJ', epsilon, E0, params{k}, ds, w(k), 10, 1e-6, [150 100]); 
    save batch.mat res params omega E H err eps
end

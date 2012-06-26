for ds = [3, 2, 1]
    for k = 1 : 4
        amp = zeros(4, 1);
        amp(k) = 1;
        [omega, E, H, err] = arm_eig('iso', amp, ds, 0.0766, epsilon, 10, 1e-6);
        res{k} = {omega, E, H, err};
    end
end
save quad_script_resd.mat res

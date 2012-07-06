function [E0] = arm_E0(res, cutout_dims)
    s = round(cutout_dims/2);
    % v = load('quad_script_resd.mat');
    for k = 1 : 4
        E = res{k}{2};
        c = round(size(E{1})/2);
        for l = 1 : 3 % cutout.
            E0{k}{l} = zeros(size(E{l}));
            E0{k}{l}(c(1)+[-s(1):s(1)], c(2)+[-s(2):s(2)], c(3)+[-s(3):s(3)]) = ...
                E{l}(c(1)+[-s(1):s(1)], c(2)+[-s(2):s(2)], c(3)+[-s(3):s(3)]);
        end

        % Find normalization factor.
        [temp, ind] = max(abs(E0{k}{2}(:)));
        norm_fact = 1 / E0{k}{2}(ind);

        for l = 1 : 3 % normalize, only keep real part.
            E0{k}{l} = real(E0{k}{l} * norm_fact);
        end
    end



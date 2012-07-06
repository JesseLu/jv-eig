function simple_movie(img, N, delay, contrast)

phase = 0;
dp = 2 * pi / N;
c = max(abs(img(:))) / contrast;
while true
    imagesc(real(img * exp(i*phase))', c * [-1 1]); axis equal tight;
    phase = phase + dp;
    pause(delay)
end


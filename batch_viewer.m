function batch_viewer(res)

subplot 111;
for k = 1 : length(res.omega)
    zc = round(size(res.E{1}{2},3)/2);
    h = imagesc(real(res.E{k}{2}(:,:,zc))'); axis equal tight;
    set(h, 'AlphaData', res.eps{k}{2}(:,:,zc)' > 5);
    res.omega{k}
    pause
end
    

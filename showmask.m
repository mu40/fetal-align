% Show 3D mask and optional points on each slice of 3D image.
function showmask(vol, mask, points, showdist)

if nargin() < 2
    points = [];
end
numpoint = size(points, 1);

reset(clf());
colormap gray
vol = stretchcon(vol);
for i = 1:size(vol, 3)
    figure(gcf());
    im = vol(:,:,i);
    if nargin() > 1
        im(mask(:,:,i)) = 1;
    end
    imagesc(im);
    for j = 1:numpoint
        pt = points(j,:);
        if abs(pt(3)-i) <= showdist
            hold on
            plot(pt(2), pt(1), 'r+', 'linewidth', 2 ,'markersize', 16);
            hold off
        end
    end
    title(sprintf('%d/%d', i, size(vol,3)));
    axis image off
    pause(0.1);
end

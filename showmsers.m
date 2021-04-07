% Show 2D masks and fitted-ellipse stats on each slice of 3D image.
function showmsers(vol, bw, snum, inplanerefdiam, subset)

if nargin() > 4
    bw = bw(subset);
    snum = snum(subset);
end
sz = size(vol);

[cen,semi,rot,~,dm] = fitellipse(bw);
ratio = semi(:,1) ./ semi(:,2);
area = prod(2*semi/inplanerefdiam, 2); % Fraction of expected size.
diam = 2*semi(:,1) / inplanerefdiam;
format = 'Slice: %2d, ratio: %.2f, area: %.2f, dice: %.2f, semi: %.2f\n';

reset(clf());
colormap gray
phi = linspace(0, 2*pi, 50);
[cx,cy] = deal(cos(phi), sin(phi));
vol = stretchcon(vol);
for i = 1:sz(3)
    figure(gcf());
    im = vol(:,:,i);
    imagesc(im);
    axis image off
    title(sprintf('%d/%d', i, sz(3)));
    drawnow();
    for j = find(snum==i)'
        tmp = im;
        tmp(bw{j}) = 1;
        imagesc(tmp);
        x = rot{j}(1,1)*semi(j,1)*cx + rot{j}(1,2)*semi(j,2)*cy + cen(j,1);
        y = rot{j}(2,1)*semi(j,1)*cx + rot{j}(2,2)*semi(j,2)*cy + cen(j,2);
        hold on
        plot(y, x, 'r', 'linewidth', 2);
        hold off
        axis('image', 'off', [0 sz(2) 0 sz(1)]);
        title(sprintf('%d/%d', i, sz(3)));
        fprintf(format, i, ratio(j), area(j), dm(j), diam(j));
        drawnow();
    end
    if ~isempty(j)
        fprintf('\n');
    end
    pause(0.1);
end

% Plot semi-major axes before/after outlier rejection.
function showfit(snum, major, zscores, upperz, coef, bound)

rejind = zscores > upperz;
[hf,hb1,hb2] = deal([]);

reset(clf());
figure(gcf());
subplot 211
hp(1) = plot(snum, major, 'bo-');
hold on
hr = plot(snum(rejind), major(rejind), 'ro');
hold off
if nargin() > 5
    hb1 = line(bound([1 1]), ylim(), 'color', 'k');
    hb2 = line(bound([2 2]), ylim(), 'color', 'k');
    lim = axis();
    xfit = lim(1):lim(2);
    hold on
    hf = plot(xfit, polyval(coef,xfit), 'g');
    hold off
    axis(lim);
end
ylabel('Semi-major axis (voxel)');

subplot 212
hp(2) = plot(snum, zscores, 'bo-');
hold on
plot(snum(rejind), zscores(rejind), 'ro');
hold off
hl = line(xlim(), [upperz upperz], 'color', 'r');
ylabel('Distance to fit (z-score)');
xlabel('Slice number');

lab = {};
if ~isempty(hr)
    lab{end+1} = 'rejected ellipses';
end
if ~isempty(hf)
    lab(end+1:end+2) = {'used for fit' 'polynomial fit'};
end
lab{end+1} = 'rejection threshold';
legend(gca(), [hr hb1 hf hl], lab, 'location', 'southeast');
set([hp hr hb1 hb2 hf hl], 'linewidth', 1.5);
drawnow();

function [clust,numclust] = meanshift(points, fwhm, gridstep, steptol, minsep)


if nargin() < 3
    gridstep = 2 * fwhm;
end
if nargin() < 4
    steptol = 0.1;
end
if nargin() < 5
    minsep = 1;
end
numdim = size(points, 2);

bounds = [min(points,[],1); max(points,[],1)];
extent = diff(bounds);
gridsize = round(extent./gridstep) + 1;
gap = extent - gridsize.*gridstep;
start = bounds(1,:) + 0.5*gap;
stop = start + gridsize.*gridstep;
clust = ndarray(start, stop, gridstep);
numclust = size(clust, 1);

sigma = fwhm / (2*sqrt(2*log(2)));
kernel = @(x,xbar)exp(-sum((x-xbar).^2,2) / (2*sigma^2));

n = 0;
while n < 100
    maxstep = 0;
    for i = 1:numclust
        weights = kernel(points, clust(i,:));
        weights = repmat(weights, [1 numdim]);
        centroid = sum(points.*weights) ./ sum(weights);
        step = sqrt(sum((centroid-clust(i,:)).^2));
        maxstep = max(step, maxstep);
        clust(i,:) = centroid;
    end
    n = n + 1;
    if maxstep < steptol
        break
    end
end

% Remove near-duplicate clusters.
i = 1;
while i <= numclust
    ind = sqrt(sum((clust(i,:)-clust).^2,2)) < minsep;
    ind(i) = 0;
    clust(ind,:) = [];
    numclust = size(clust, 1);
    i = i + 1;
end

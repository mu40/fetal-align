function [clusters,numclust] = meanshift(points, diam, steptol, minsep)

if nargin() < 3
    steptol = 0.1;
end
if nargin() < 4
    minsep = 1;
end
numdim = size(points, 2);

bounds = [min(points,[],1); max(points,[],1)];
extent = diff(bounds);
gridsize = round(extent./diam) + 1;
gap = 0.5 * (extent-gridsize.*diam);
start = bounds(1,:) + 0.5*gap;
stop = start + gridsize.*diam;
clusters = ndarray(start, stop, diam);
numclust = size(clusters, 1);

fwhm = diam/2;
sigma = fwhm / (2*sqrt(2*log(2)));
kernel = @(x,xbar)exp(-sum((x-xbar).^2,2) / (2*sigma^2));

n = 0;
while n < 100
    maxstep = 0;
    for i = 1:numclust
        weights = kernel(points, clusters(i,:));
        weights = repmat(weights, [1 numdim]);
        centroid = sum(points.*weights) ./ sum(weights);
        step = sqrt(sum((centroid-clusters(i,:)).^2));
        maxstep = max(step, maxstep);
        clusters(i,:) = centroid;
    end
    n = n + 1;
    if maxstep < steptol
        break
    end
end

% Remove near-duplicate clusters.
i = 1;
while i <= numclust
    ind = sqrt(sum((clusters(i,:)-clusters).^2,2)) < minsep;
    ind(i) = 0;
    clusters(ind,:) = [];
    numclust = size(clusters, 1);
    i = i + 1;
end

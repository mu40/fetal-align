% Fit line to ND point cloud: return centroid, direction and point distances.
% http://www.ctralie.com/Teaching/COMPSCI290/Lectures/10_PCA/slides.pdf
function [centroid, dir, dist] = fitline(points)

centroid = mean(points, 1);
points = points - centroid;

covmat = points' * points;
[evec,eval] = eig(covmat, 'vector');

[~,order] = sort(abs(eval), 'descend');
evec = evec(:,order);
dir = evec(:,1)';

p = points - centroid;
dist = real(sqrt(sum(p.*p,2)-(sum(p.*dir,2)).^2));

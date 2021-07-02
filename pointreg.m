% Find transform from point cloud X to Y: Umeyama et al., 1991. The sets must
% be of equal size and the points must correspond across sets X and Y.'''
function [rot,tra,sca,err] = pointreg(setx, sety, scaling)

if nargin() < 3
	scaling = 0;
end

numdim = size(setx, 1);
numpts = size(setx, 2);
meanx = mean(setx, 2);
meany = mean(sety, 2);
varx = sum(sum((setx-meanx).^2, 1)) / numpts;
vary = sum(sum((sety-meany).^2, 1)) / numpts;
covmat = (sety - meany) * (setx - meanx)' / numpts;
[u,d,v] = svd(covmat);

s = eye(numdim);
if det(covmat) < 0
	s(end) = -1;
end
rots = s;
if rank(covmat) == numdim - 1
	rots = eye(numdim);
	if round(det(u)*det(v)) < 0
		rots(end) = -1;
	end
end
rot = u * rots * v';

if scaling
	sca = trace(d*s) / varx;
else
	sca = 1;
end
tra = meany - sca * rot * meanx;
err = sqrt(abs(vary - trace(d*s)^2 / varx));

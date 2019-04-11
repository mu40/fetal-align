% Generate coordinates of ellipsoid surface for use with MATLAB's surf()
% function. By specifying a permution of [1 2 3] as order, the poles of the
% resulting mesh can be placed somewhere else than at extremal z (3).
function [x,y,z] = genellipsoid(center, radius, rotmat, degstep, order)

if nargin() > 4 && ~all(ismember([1 2 3], order))
    error('Order isn''t a permutation of [1 2 3].');
end
if nargin() < 4
    degstep = 10;
end
if numel(radius) == 1
    radius = repmat(radius, [1 3]);
end
if isempty(rotmat)
    rotmat = eye(3);
end

step = pi/180 * degstep;
theta = 0:step:pi;
phi = 0:step:2*pi;
numtheta = numel(theta);
numphi = numel(phi);

[x,y,z] = deal(zeros(numtheta, numphi));
for i = 1:numtheta
    x(i,:) = sin(theta(i)) .* cos(phi);
    y(i,:) = sin(theta(i)) .* sin(phi);
    z(i,:) = repmat(cos(theta(i)), [1 numphi]);
end

if nargin() > 4
    tmp = cat(3, x, y, z);
    x = tmp(:,:,order(1));
    y = tmp(:,:,order(2));
    z = tmp(:,:,order(3));
end

xyz = [x(:) y(:) z(:)];
xyz = xyz .* radius;
if nargin() > 2
    xyz = (rotmat * xyz')';
end
xyz = xyz + center;

x = reshape(xyz(:,1), size(x));
y = reshape(xyz(:,2), size(y));
z = reshape(xyz(:,3), size(z));

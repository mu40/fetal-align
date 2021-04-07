% Fit an ellipse/ellipsoid to 2D/3D image blobs.
function [bcen,semilen,rotmat,ellmask,dm] = ...
    fitellipse(mask, scalemask)

if nargin() < 2
    scalemask = 1;
end

iscellinput = iscell(mask);
if ~iscellinput
    mask = {mask};
end
sz = single(size(mask{1}));
numdim = numel(sz);
if numel(sz) == 2
    sz(3) = 1;
end
[x,y,z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));

nummask = numel(mask);
bcen = zeros(nummask, numdim);
semilen = zeros(nummask, numdim);
rotmat = cell(nummask, 1);
ellmask = cell(nummask, 1);
dm = zeros(nummask, 1);
for i = 1:nummask
    % Centre of mass and barycentric coordinates.
    bw = mask{i};
    ind = find(bw);
    m000 = numel(ind);
    cen = [sum(x(ind)) sum(y(ind)) sum(z(ind))] / m000;
    xb = x(ind) - cen(1);
    yb = y(ind) - cen(2);
    zb = z(ind) - cen(3);
    
    % Second moments - weighted by intensity, which is 1.
    m200 = sum(xb.^2);
    m020 = sum(yb.^2);
    m002 = sum(zb.^2);
    m110 = sum(xb.*yb);
    m011 = sum(yb.*zb);
    m101 = sum(xb.*zb);
    mi = [m200 m110 m101; m110 m020 m011; m101 m011 m002] / m000;
    
    % Semi axes lengths and rotation matrix from eigendecomposition.
    [evec,eval] = eig(mi, 'vector');
    [eval,order] = sort(eval, 'descend');
    rot = evec(:,order);
    semi = 2 * real(sqrt(eval)); % Tiny imaginary part of 2D structure in 3D.
    
    % Mask of ellipse.
    halflen = semi(1);
    lower = max(1, floor(cen-halflen));
    upper = min(sz, ceil(cen+halflen));
    indx = lower(1):upper(1);
    indy = lower(2):upper(2);
    indz = lower(3):upper(3);
    xb = x(indx,indy,indz) - cen(1);
    yb = y(indx,indy,indz) - cen(2);
    zb = z(indx,indy,indz) - cen(3);
    xr = rot(1,1)*xb + rot(2,1)*yb + rot(3,1)*zb;
    yr = rot(1,2)*xb + rot(2,2)*yb + rot(3,2)*zb;
    zr = rot(1,3)*xb + rot(2,3)*yb + rot(3,3)*zb;
    ell = false(sz);
    scale = semi .* scalemask;
    scale(scale == 0) = 1; % Avoid 2D division by zero.
    isell = (xr./scale(1)).^2 + (yr./scale(2)).^2 + (zr./scale(3)).^2 < 1;
    ell(indx,indy,indz) = isell;
    
    if sz(3) == 1
        cen = cen(1:2);
        semi = semi(1:2);
        rot = rot(1:2,1:2);
    end
    rot = sign(det(rot)) * rot; % Proper rotation matrix.
    
    bcen(i,:) = cen;
    semilen(i,:) = semi;
    rotmat{i} = rot;
    ellmask{i} = ell;
    dm(i) = dice(bw, ell);
end

if ~iscellinput
    rotmat = rotmat{1};
    ellmask = ellmask{1};
end

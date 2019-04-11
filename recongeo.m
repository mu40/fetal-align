% Reconstruct fetal cordinate system and return tranfrom from it to RAS space.
function fetaltoras = recongeo(voxtoras, bcenvox, ecenvox, bmask)

% Matrix that rotates x-y-z onto major-middle-minor axes.
[~,semivox,rotvox] = fitellipse(bmask);
voxsize = sqrt(sum(voxtoras(1:3,1:3).^2));
evecmm = diag(voxsize) * rotvox * diag(semivox);
semimm = sqrt(sum(evecmm.^2));
[~,order] = sort(semimm, 'desc');
rotvox = rotvox(:,order);
rotvox = sign(det(rotvox)) * rotvox; % Force proper rotation matrix.

% Convert to world space.
bcenmm = voxtoras * [bcenvox 1]';
ecenmm = voxtoras * [ecenvox'; 1 1];
bcenmm = bcenmm(1:3);
ecenmm = ecenmm(1:3,:);

% Anatomical axes: derive HF from LR and AP, re-estimate AP (first LR/AP not
% exactly perpendicular), re-estimate LR - to obtain proper rotation matrix.
lr = ecenmm(:,2) - ecenmm(:,1);
pa = 0.5 * (ecenmm(:,2)-bcenmm + ecenmm(:,1)-bcenmm);
lr = lr / norm(lr);
pa = pa / norm(pa);
fh = cross(lr, pa); % RAS space.
fh = fh / norm(fh);
lr = cross(pa, fh);

% Orientation of HF axis: plane LR/major ellipsoid axis divides brain into top
% and bottom halves. Normal vector of that plane points to half space with the
% eyes if HF oriented correctly.
major = voxtoras(1:3,1:3) * rotvox * [1 0 0]';
major = major - lr*dot(major, lr);
major = major / norm(major);
normal = cross(major, lr);
normal = normal / norm(normal);
if dot(normal, pa) < 0
    normal = -1 * normal; % Now pointing to side of eyes.
end
if dot(major, pa) < 0
    major = -1 * major; %#ok % Now pointing forward.
end
if dot(normal, fh) > 0 % HF pointing towards side of eyes, incorrecly.
    fh = -fh;
    lr = -lr;
end

% Transform from brain to world.
fetaltoras = eye(4);
fetaltoras(1:3,1:3) = [lr pa fh];
incline = rotmat(30 .* pi/180); % Rotate from eyes around LR axis.
shift = [eye(4,3) [bcenmm; 1]];
fetaltoras = shift * fetaltoras * incline;

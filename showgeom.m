% Display geometry of brain and eyes in scaled voxel space.
function showgeom(brainmask, ecenvox, voxsize, odiam)

% Scaled voxel space.
[bcenvox,bsemivox,rotvox] = fitellipse(brainmask);
evecmm = diag(voxsize) * rotvox * diag(bsemivox);
bsemimm = sqrt(sum(evecmm.^2));
[bsemimm,order] = sort(bsemimm, 'desc');
rotmm = rotvox(:,order);
rotmm = sign(det(rotmm)) * rotmm; % Force proper rotation matrix.
bcenmm = bcenvox .* voxsize;
ecenmm = ecenvox .* voxsize;
numeye = size(ecenmm, 1);

% Geometry.
[bx,by,bz] = genellipsoid(bcenmm, bsemimm, rotmm, 10, [3 2 1]);
[ex,ey,ez] = genellipsoid([0 0 0], odiam/2, rotmm, 30, [3 2 1]);
major = 1.2*bsemimm(1) .* rotmm * [1 0 0]';
major = bcenmm' + [-1 1] .* major;
[lx,ly,lz] = deal(major(1,:), major(2,:), major(3,:));

reset(clf());
figure(gcf());
hb = surf(bx, by, bz);
hold on
for i = 1:numeye
    he(i) = surf(ecenmm(i,1)+ex, ecenmm(i,2)+ey, ecenmm(i,3)+ez); %#ok
end
hl = line(lx, ly, lz);
hold off
axis equal
legend([hb hl he(1)], {'brain', 'major axis' 'eyes'});
xlabel('Phase (mm)');
ylabel('Read (mm)');
zlabel('Slice (mm)');
set(hb, 'facecolor' , 'k', 'edgecolor', 'k', 'edgealpha', 0.2);
set(he(1:min(2,end)), 'facecolor', 'r', 'edgecolor', 'r');
set(he(3:end), 'facecolor', 'b', 'edgecolor', 'b');
set([hb he], 'facealpha', 0.1);
set(hl, 'linewidth', 2, 'color', 'k');
set(gca, 'projection', 'perspective', 'cameratarget', bcenmm);

% Movement.
theta = 70;
phi = -37.5 + (0:1:180);
ang = get(gca, 'cameraviewangle');
target = get(gca, 'cameratarget');
camera = get(gca, 'cameraposition');
radius = norm(camera - target);
for i = 1:numel(phi)
    figure(gcf());
    sx = sind(theta) .* cosd(phi(i));
    sy = sind(theta) .* sind(phi(i));
    sz = cosd(theta);
    pos = radius .* [sx sy sz] + target;
    set(gca, 'cameraposition', pos, 'cameraviewangle', ang);
    pause(0.01);
end

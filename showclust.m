function showclust(points, clusters, diam)

numclust = size(clusters, 1);
[sx,sy,sz] = genellipsoid([0 0 0], diam/2, eye(3), 10);


reset(clf());
figure(gcf());
hm = plot3(points(:,1), points(:,2), points(:,3), 'r+');
hold on
hp = plot3(clusters(:,1), clusters(:,2), clusters(:,3), 'b+');
for i = 1:numclust
    hs(i) = surf(sx+clusters(i,1), sy+clusters(i,2), sz+clusters(i,3)); %#ok
end
hold off
grid on
axis equal
xlabel('Phase (mm)');
ylabel('Read (mm)');
zlabel('Slice (mm)');
legend('MSER centers', 'cluster centers', 'd=OFD spheres');
set(hs, 'facecolor' , 'k', 'edgecolor', 'k', 'edgealpha', 0.2, ...
    'facealpha', 0.1);
set([hp hm], 'linewidth', 2);
set(gca, 'projection', 'perspective');

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

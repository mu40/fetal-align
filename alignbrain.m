% Reslice brain such that it's aligned with the anatomical planes, based on
% the coordinates of its centre and those of the eyes (voxel coordinates).
function out = alignbrain(mri, fetaltoras, outfile)

if nargin() < 3
    outfile = [];
end

if ischar(mri)
    mri = MRIread(mri); % FreeSurfer.
    mri.vol = permute(mri.vol, [2 1 3 4]);
end
voxtoras = mri.vox2ras1;
dat = mri.vol(:,:,:,1);
dat = dat ./ gaussblur(dat, 0.33*size(dat));

% Resample in anatomical coordinates.
fov = -70:1:70;
[xf,yf,zf] = ndgrid(fov, fov, fov);
fetaltovox = voxtoras \ fetaltoras;
[xi,yi,zi] = coords(fetaltovox, xf, yf, zf);
res = interpn(dat, xi, yi, zi, 'linear', 0);

% Plot central slices.
reset(clf());
figure(gcf());
colormap gray

hs(1) = subplot(2,2,1);
im = squeeze(res(:,round(end/2),:));
im = flip(im', 1);
imagesc(stretchcon(im));
title('Coronal');

hs(2) = subplot(2,2,2);
im = squeeze(res(round(end/2),:,:));
im = flip(flip(im', 1),2);
imagesc(stretchcon(im));
title('Sagittal');

hs(3) = subplot(2,2,3);
im = squeeze(res(:,:,round(end/2)));
im = flip(im', 1);
imagesc(stretchcon(im));
title('Axial');

hl = [];
for i = 1:3
    hl(end+1) = line(hs(i), xlim(), [1 1].*mean(ylim())); %#ok<*AGROW>
    hl(end+1) = line(hs(i), [1 1].*mean(xlim()), ylim());
end
axis(hs, 'image', 'off');
set(hl, 'color', 'y', 'linewidth', 1.5);
drawnow();

out = struct();
out.vol = permute(res, [2 1 3 4]);
voxsize = [xf(2,1,1)-xf(1,1,1) yf(1,2,1)-yf(1,1,1) zf(1,1,2)-zf(1,1,1)];
out.volres = voxsize;
scaling = diag([voxsize 1]);
subone = [eye(4,3) [-1 -1 -1 1]'];
voxtofetal = [eye(4,3) [xf(1) yf(1) zf(1) 1]'] * scaling * subone;
out.vox2ras1 = fetaltoras * voxtofetal;
out.vox2ras0 = vox2ras_1to0(out.vox2ras1);
if ~isempty(outfile)
    MRIwrite(out, outfile);
end

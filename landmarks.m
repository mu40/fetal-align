% Estimate center of fetal brain and largest slices from 2D-EPI scout.
function [bcenout,ecenout,outmask] = landmarks(mri, ga, doplot, par)

if nargin() < 3
    doplot = [];
end
if nargin() < 4
    par = [];
end

if ischar(mri)
    mri = MRIread(mri); % FreeSurfer.
    mri.vol = permute(mri.vol, [2 1 3 4]);
end
dat = single(mri.vol(:,:,:,1));
vsz = mri.volres;
dim = mri.volsize;

crop = 0;
if any(dim > 128)
    crop = 0.1;
end
cropshift = floor([crop*dim(1:2) 0]);
[low,upp] = deal(1+cropshift, dim-cropshift);
dat = dat(low(1):upp(1), low(2):upp(2), low(3):upp(3));
dim = size(dat);

% Pre-processing.
par = setdefault(par, 'nufwhm', 20, 'dogfwhm', 5, 'dogratio', 2);
[nuc, pre, flt] = deal(zeros(dim, 'single'));
fwhm = par.dogfwhm ./ vsz(1:2);
for i = 1:dim(3)
    im = dat(:,:,i);
    im = im ./ gaussblur(im, par.nufwhm./vsz(1:2));
    nuc(:,:,i) = im;
    tmp = gaussblur(im, par.dogratio*fwhm) - gaussblur(im, fwhm);
    tmp(tmp<0) = 0;
    tmp = 1 - stretchcon(tmp);
    flt(:,:,i) = tmp;
    im = im .* tmp;
    im = stretchcon(im, [10 99]);
    pre(:,:,i) = im;
end
[ofd,bpd,odiam,odist] = anatomy(ga);

%% Stage 1.

% MSERs and filtering.
par = setdefault(par, 'maxratio1', 1.5, 'maxarea1', 1.1, 'maxsemi1', 1.1, ...
    'mindice1', 0.7, 'minmser1', 0.2, 'maxmser1', 1, 'delta', 5);
refarea = pi * bpd/2 * ofd/2 / prod(vsz(1:2));
[bw, numbw, snum] = slcmser(pre, par.delta, refarea*par.minmser1, ...
    refarea*par.maxmser1);
[cen,semi,rot,~,dice] = fitellipse(bw); %#ok<ASGLU>
ratio = semi(:,1) ./ semi(:,2);
area = pi * prod(semi, 2);
keep = true(numbw, 1);
keep = keep & ratio < ofd/bpd * par.maxratio1;
keep = keep & area < refarea * par.maxarea1;
keep = keep & semi(:,1)*vsz(1) < ofd/2 * par.maxsemi1;
keep = keep & dice > par.mindice1;

% Mean-shift clustering: remove near-duplicate points first.
par = setdefault(par, 'clustfwhm1', 1, 'gridstep1', 1, 'minsep1', 1, ...
    'steptol1', 0.1);
points = [cen(keep,:) snum(keep)] .* vsz;
i = 1;
while i <= size(points,1)
    ind = sqrt(sum((points(i,:)-points).^2,2)) < par.minsep1; % In mm.
    ind(i) = 0;
    points(ind,:) = [];
    i = i + 1;
end
[clusters,numclust] = meanshift(points, par.clustfwhm1*ofd/2, ...
    par.gridstep1*ofd, par.steptol1, par.minsep1);
if isfield(doplot, 'clust1') && doplot.clust1
    wait('Brain localization: mean-shift clustering of MSER centers', ...
        'before cluster selection');
    showclust(points, clusters, par.gridstep1*ofd);
end

% Cluster selection.
par = setdefault(par, 'clustrad1', 1);
radius = ofd/2 * par.clustrad1;
xyz = ndarray([1 1 1], dim) .* vsz;
ind1d = reshape(1:prod(dim), dim);
points = [cen snum] .* vsz;
[numvoxin, numvoxout] = deal(zeros(numclust, 1));
for i = 1:numclust
    centroid = clusters(i,:);
    dist = sqrt(sum((centroid-points).^2, 2));
    ind = find(keep & dist<radius);
    volmask = false(dim);
    for j = ind'
        volmask(:,:,snum(j)) = volmask(:,:,snum(j)) | bw{j};
    end
    low = floor((centroid-radius) ./ vsz);
    upp = ceil((centroid+radius) ./ vsz);
    low = max(1, low);
    upp = min(dim, upp);
    subset = ind1d(low(1):upp(1),low(2):upp(2),low(3):upp(3));
    sphere = false(dim);
    insphere = sqrt(sum((xyz(subset,:)-centroid).^2, 2)) < radius;
    sphere(subset) = insphere;
    numvoxin(i) = nnz(sphere & volmask);
    numvoxout(i) = nnz(~sphere & volmask);
end
score = numvoxin - numvoxout;
[~,order] = sort(score, 'descend');
centroid = clusters(order(1),:);
dist = sqrt(sum((centroid-points(keep,:)).^2, 2));
keeppreclust = keep; %#ok
keep(keep) = dist < radius;
if isfield(doplot, 'mser1') && doplot.mser1
    wait('Brain localization: retained MSERs on each slice before', ...
        'combination into a preliminary brain mask');
    showmsers(pre, bw, snum, ofd/vsz(1), keep);
end

% Create mask and find barycenter.
rawmask = false(dim);
for i = find(keep)'
    rawmask(:,:,snum(i)) = rawmask(:,:,snum(i)) | bw{i};
end
for i = 1:dim(3)
    rawmask(:,:,i) = imfill(rawmask(:,:,i), 4, 'holes');
end
rawmask = closegaps(rawmask);
rawmask = cutprot(rawmask);
rawloc = fitellipse(rawmask);
if isfield(doplot, 'mask1') && doplot.mask1
    wait('Brain localization: slices of the preliminary brain mask');
    showmask(pre, rawmask, rawloc, ofd/vsz(3)/2);
end

%% Stage 2.

% Image cropping.
halflen = ofd ./ vsz;
low = max(1, floor(rawloc-halflen));
upp = min(dim, floor(rawloc+halflen));
boxshift = low - 1;
boxdat = pre(low(1):upp(1), low(2):upp(2), low(3):upp(3));
boxdim = size(boxdat);
rawloc = rawloc - boxshift;

% MSER detection and filtering.
par = setdefault(par, 'maxratio2', 1.5, 'maxarea2', 1.1, 'maxsemi2', 1.2, ...
    'mindice2', 0.7, 'maxdist2', 1.1, 'minmser2', 0.05, ...
    'maxmser2', 1);
refarea = pi * bpd/2 * ofd/2 / prod(vsz(1:2));
[bw, numbw, snum] = slcmser(boxdat, par.delta, refarea*par.minmser2, ...
    refarea*par.maxmser2);
[cen,semi,~,~,dice] = fitellipse(bw);
ratio = semi(:,1) ./ semi(:,2);
area = pi * prod(semi, 2);
dist = sqrt(sum((vsz.*([cen snum]-rawloc)).^2, 2));
keep = true([numbw 1]);
keep = keep & ratio < ofd/bpd * par.maxratio2;
keep = keep & area < refarea * par.maxarea2;
keep = keep & semi(:,1)*vsz(1) < ofd/2 * par.maxsemi2;
keep = keep & dice > par.mindice2;
keep = keep & dist < ofd/2 * par.maxdist2;
if isfield(doplot, 'mser2') && doplot.mser2
    wait('Brain-mask creation: pre-filtered MSERs on each slice before', ...
        'addition to the preliminary brain mask');
    showmsers(boxdat, bw, snum, ofd/vsz(1), keep);
end

% Brain mask update.
volmask = rawmask(low(1):upp(1), low(2):upp(2), low(3):upp(3));
for i = find(keep)'
    volmask(:,:,snum(i)) = volmask(:,:,snum(i)) | bw{i};
end
% for i = 1:boxdim(3)
%     volmask(:,:,i) = imfill(volmask(:,:,i), 4, 'holes');
% end
volmask = closegaps(volmask);
boxloc = fitellipse(volmask);

% Slice-wise ellipses.
masksize = sum(sum(volmask));
slcind = find(masksize);
slcmasks = num2cell(volmask(:,:,slcind), [1 2]);
[scen,ssemi] = fitellipse(slcmasks);
scen = [scen slcind];

% Fit polynomial to average axes.
[scenprepoly, ssemiprepoly, boxlocprepoly] = deal(scen, ssemi, boxloc); %#ok
par = setdefault(par, 'polywithin2', 0.5, 'zpoly2', 1.4);
polywithin = ofd/vsz(3)/2 * par.polywithin2;
ind = abs(scen(:,3)-boxloc(3)) < polywithin;
if nnz(ind) > 2 % Need >=3 points for quadratic fit.
    xdat = scen(ind,3);
    ydat = ssemi(ind,1);
    coef = polyfit(xdat, ydat, 2);
    dist = ssemi(:,1) - polyval(coef, scen(:,3));
    zscorepoly = (dist-mean(dist(ind))) / std(dist);
    if coef(1) < 0 % Want bad weather.
        ind = zscorepoly < par.zpoly2;
        volmask(:,:,scen(~ind,3)) = 0;
        scen = scen(ind,:); %#ok
        ssemi = ssemi(ind,:); %#ok
    end
    if isfield(doplot, 'poly2') && doplot.poly2
        wait('Brain-mask creation: rejection of brain-mask slices', ...
            'based on the deviation between their major-axis length', ...
            'and a quadratic fit across slices');
        showfit(scenprepoly(:,3), ssemiprepoly(:,1), zscorepoly, ...
            par.zpoly2, coef, boxloc(3)+[-1 1]*polywithin);
    end
end

% Aggressive final gap filling: included anywhere before and after.
volmask = cutprot(volmask);
for i = 2:boxdim(3)-1
    before = sum(volmask(:,:,1:i-1), 3);
    after = sum(volmask(:,:,i+1:end), 3);
    volmask(:,:,i) = volmask(:,:,i) | (before & after);
end

% Outputs in full image space.
[x,y,z] = ndgrid(1:boxdim(1), 1:boxdim(2), 1:boxdim(3));
ind = find(volmask);
boxloc = [sum(x(ind)) sum(y(ind)) sum(z(ind))] / numel(ind);
bcenvox = boxloc + boxshift;
brainmask = false(dim);
brainmask(low(1):upp(1), low(2):upp(2), low(3):upp(3)) = volmask;
if isfield(doplot, 'mask2') && doplot.mask2
    wait('Brain-mask creation: slices of the final brain mask');
    showmask(boxdat, volmask, boxloc, ofd/vsz(3)/2);
end

%% Stage 3.

% Image cropping and interpolation.
par = setdefault(par, 'voxsize3', 1.49);
boxvsz = [par.voxsize3 par.voxsize3 vsz(3)];
lenvox = floor(sqrt(2) * ofd ./ boxvsz); % Side length in new voxels.
scaledown = boxvsz ./ vsz;
shift = bcenvox - (lenvox+1)/2 .* scaledown; % In old voxels.
shift(3) = round(shift(3)); % Don't average slices (motion).
shift = max(shift, [0 0 0]); % Don't start before voxel 1.
shift = [eye(4,3) [shift 1]'];
scaledown = diag([scaledown 1]);
boxtofull = shift * scaledown;
maxlenvox = floor(boxtofull \ [dim 1]');
maxlenvox = maxlenvox(1:3)';
lenvox = min(lenvox, maxlenvox); % Don't go beyond last voxel.
[xb,yb,zb] = ndgrid(1:lenvox(1), 1:lenvox(2), 1:lenvox(3));
[xb,yb,zb] = coords(boxtofull, xb, yb, zb);
boxdat = interpn(pre, xb, yb, zb, 'linear', 0);
boxdim = size(boxdat);
boxloc = boxtofull \ [bcenvox 1]';
boxloc = boxloc(1:3)';

% MSER detection and filtering.
par = setdefault(par, 'maxratio3', 1.5, 'maxarea3', 1.2, 'maxsemi3', 1.3, ...
    'mindice3', 0.9, 'mindist3', 0.6, 'maxdist3', 1.3, 'mindrop3', 0.5, ...
    'ribbon3', 1.5, 'minmser3', 0.5, 'maxmser3', 1.1);
refarea = pi * (odiam/2)^2 / prod(boxvsz(1:2)); % In voxels.
[bw, numbw, snum] = slcmser(boxdat, par.delta, refarea*par.minmser3, ...
    refarea*par.maxmser3);
[cen,semi,rot,ell,dice] = fitellipse(bw); %#ok
ratio = semi(:,1) ./ semi(:,2);
area = pi * prod(semi, 2);
[~,~,~,ellbig] = fitellipse(bw, par.ribbon3);
contrast = zeros(numbw, 1);
for i = 1:numbw
    im = boxdat(:,:,snum(i));
    inribbon = find(ellbig{i} & ~ell{i});
    contrast(i) = median(im(inribbon(:))) / median(im(ell{i}(:)));
end
dist = sqrt(sum((boxvsz .* ([cen snum]-boxloc)).^2, 2));
keep = true(numbw, 1);
keep = keep & ratio < par.maxratio3;
keep = keep & area < refarea * par.maxarea3;
keep = keep & semi(:,1).*boxvsz(1) < odiam/2 * par.maxsemi3;
keep = keep & dice > par.mindice3;
keep = keep & contrast < par.mindrop3;
keep = keep & dist > ofd/2 * par.mindist3;
keep = keep & dist < ofd/2 * par.maxdist3;
if isfield(doplot, 'mser3') && doplot.mser3
    wait('Eye detection: pre-filtered MSERs on each slice before 3D', ...
        'clustering');
    showmsers(boxdat, bw, snum, ofd/boxvsz(1), keep);
end

% Clusters of points in 3D.
par = setdefault(par, 'group3', 0.8);
ind = find(keep);
cluster = arrayfun(@(p)p, ind, 'uniformoutput', 0);
numclust = numel(cluster);
for i = 1:numclust
    cen1 = [cen(ind(i),:) snum(ind(i))] .* boxvsz;
    for j = 1:numclust
        cen2 = [cen(ind(j),:) snum(ind(j))] .* boxvsz;
        if i ~= j && sqrt(sum((cen1-cen2).^2)) < odiam*par.group3
            cluster{i} = [cluster{i} ind(j)];
        end
    end
end
weight = rand(1, numclust);
hash = cellfun(@(x)sum(weight(1:numel(x)).*sort(x,'ascend')), cluster);
[~,ind] = unique(hash);
cluster = cluster(ind);
numclust = numel(ind);

% Cluster: center, number of slices, quality.
ccen = zeros([numclust 3]);
numslc = zeros([numclust 1]);
[x,y,z] = ndgrid(1:boxdim(1), 1:boxdim(2), 1:boxdim(3));
for i = 1:numclust
    clustind = cluster{i};
    mask = false(boxdim);
    onslice = false([1 boxdim(3)]);
    for j = clustind
        slcind = snum(j);
        mask(:,:,slcind) = mask(:,:,slcind) | bw{j};
        onslice(slcind) = 1;
    end
    numslc(i) = sum(onslice);
    ismask = find(mask);
    ccen(i,:) = sum([x(ismask) y(ismask) z(ismask)]) / numel(ismask);
end

% Distance to brain center.
locmm = boxloc .* boxvsz;
ccenmm = ccen .* boxvsz;
bdismm = sqrt(sum((ccenmm-locmm).^2, 2));

% Distance between eyes.
x = ccenmm(:,1);
y = ccenmm(:,2);
z = ccenmm(:,3);
edismm = sqrt((x-x').^2 + (y-y').^2 + (z-z').^2);

% Centroid pair scores.
par = setdefault(par, 'lambda3', 4, 'theta3', 1, 'int3', 1);
numpair = numclust*(numclust-1)/2;
score = zeros([numpair 1]);
pair = zeros([numpair 2]);
n = 1;
f = @(x)abs(x);
for i = 1:numclust
    for j = i+1:numclust
        tmp = 0;
        tmp = tmp + 1.0 * f(mean(bdismm([i j]))/(ofd/2) - 1);
        tmp = tmp + par.lambda3 * f((bdismm(i)-bdismm(j))/mean(bdismm([i j])));
        tmp = tmp + par.int3 * f(edismm(i,j)/odist - 1);
        tmp = tmp + par.theta3 * f(mean(numslc([i j]))/max(numslc) - 1);
        score(n) = tmp;
        pair(n,:) = [i j];
        n = n + 1;
    end
end
[score,order] = sort(score, 'ascend'); %#ok<ASGLU>
pair = pair(order,:);
finccen = ccen(pair(1,:),:);
finccenmm = ccenmm(pair(1,:),:);
if ~isempty(doplot)
    wait('Eye detection: geometric properties after cluster selection');
    fprintf('|E1-B|/(0.5*OFD): %.2f\n', bdismm(pair(1,1)) / (0.5*ofd));
    fprintf('|E2-B|/(0.5*OFD): %.2f\n', bdismm(pair(1,2)) / (0.5*ofd));
    fprintf('|E1-E2|/ODIST:    %.2f\n', edismm(pair(1,1),pair(1,2)) / odist);
    fprintf('Cluster slices:   %d %d\n', numslc(pair(1,:)));
    fprintf('\n');
end

eyemask = false(boxdim);
dis1mm = sqrt(sum((finccenmm(1,:)-[cen snum].*boxvsz).^2, 2));
dis2mm = sqrt(sum((finccenmm(2,:)-[cen snum].*boxvsz).^2, 2));
for i = 1:boxdim(3)
    ind = keep & snum==i & (dis1mm<odiam/2 | dis2mm<odiam/2);
    for j = find(ind)'
        eyemask(:,:,i) = eyemask(:,:,i) | bw{j};
    end
end
if isfield(doplot, 'mask3') && doplot.mask3
    wait('Eye detection: slices of the eye masks after cluster selection');
    showmask(boxdat, eyemask, finccen, odiam/boxvsz(3)/2);
end

% Center of eyes in original voxel space.
ecenvox = boxtofull * [finccen'; 1 1];
ecenvox = ecenvox(1:3,:)';
if isfield(doplot, 'geom3') && doplot.geom3
    wait('Eye detection: spatial configuration of the detected landmarks');
    showgeom(brainmask, ecenvox, vsz, odiam);
end

%% Housekeeping.

% Account for cropping.
bcenout = bcenvox + cropshift;
ecenout = ecenvox + cropshift;
outmask = false(mri.volsize);
xsub = 1+cropshift(1) : mri.volsize(1)-cropshift(1);
ysub = 1+cropshift(2) : mri.volsize(2)-cropshift(2);
zsub = 1+cropshift(3) : mri.volsize(3)-cropshift(3);
outmask(xsub,ysub,zsub) = brainmask;

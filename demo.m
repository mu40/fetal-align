% Demo visualizing fetal-brain localization and orientation estimation.

% Plot stages: 1 brain localization, 2 masking, 3 eye detection (0 is off).
doplot = struct();
doplot.clust1 = 1;
doplot.line1 = 1;
doplot.mser1 = 1;
doplot.mask1 = 1;
doplot.mser2 = 1;
doplot.line2 = 1;
doplot.poly2 = 1;
doplot.mask2 = 1;
doplot.mser3 = 1;
doplot.mask3 = 1;
doplot.geom3 = 1;

% Load data paths and gestational age.
addpath freesurfer
[data,ga] = loaddata();
wait('This demo showcases the different stages of fetal-brain', ...
    'localization and orientation estimation. You will be prompted for', ...
    'input before each step is run.');

% Import image.
imnum = 1;
mri = MRIread(data{imnum}); % FreeSurfer.
mri.vol = permute(mri.vol, [2 1 3 4]);
wait('Slices of input image', ['"' data{imnum} '"']);
showmask(stretchcon(mri.vol));

% Detect landmarks and estimate orientation.
[bvox,evox,bmask] = landmarks(mri, ga(imnum), doplot);
fettoras = estorient(mri.vox2ras1, bvox, evox, bmask);

% Resample in anatomical space and save for inspection. 
outfile = 'out.nii';
wait('Resampling in anatomical space');
alignbrain(mri, fettoras, outfile);
fprintf('Done\n');

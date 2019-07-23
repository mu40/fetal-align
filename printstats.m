% Deviation between automatically and manually determined landmarks.

addpath freesurfer
[data,ga,b_manvox,e1manvox,e2manvox] = loaddata();

numdat = numel(data);
successful = setdiff(1:numdat, [34 41]); % Exclude failures.
numdat = numel(successful);
[magrot,magtra,mag_e1,mag_e2] = deal(zeros(numdat, 1));
[xtra,ytra,ztra,xrot,yrot,zrot] = deal(zeros(numdat, 1));
for i = 1:numdat
    fprintf('Processing dataset %d/%d\n', i, numdat);
    ind = successful(i);
    mri = MRIread(data{ind}); % FreeSurfer.
    mri.vol = permute(mri.vol, [2 1 3 4]);
    worldmat = mri.vox2ras1;
    [b_aut,e_aut,bmask] = landmarks(mri, ga(ind));
    % fettoras = estorient(worldmat, b_aut, e_aut, bmask);
	% alignbrain(mri, fettoras);
    
    % Transform to RAS.
    b_aut = worldmat * [b_aut 1]';
    e1aut = worldmat * [e_aut(1,:) 1]';
    e2aut = worldmat * [e_aut(2,:) 1]';
    b_man = worldmat * [b_manvox{ind} 1]';
    e1man = worldmat * [e1manvox{ind} 1]';
    e2man = worldmat * [e2manvox{ind} 1]';
    [b_aut,e1aut,e2aut] = deal(b_aut(1:3), e1aut(1:3), e2aut(1:3));
    [b_man,e1man,e2man] = deal(b_man(1:3), e1man(1:3), e2man(1:3));
    
     % Deviation for eyes. Labelling mismatch?
    if norm(e1aut-e2man) < norm(e1aut-e1man)
        tmp = e1aut;
        e1aut = e2aut;
        e2aut = tmp;
    end
    mag_e1(i) = norm(e1aut - e1man);
    mag_e2(i) = norm(e2aut - e2man);
    
    % Deviation for brain.
    xtra(i) = b_aut(1) - b_man(1);
    ytra(i) = b_aut(2) - b_man(2);
    ztra(i) = b_aut(3) - b_man(3);
    magtra(i) = sqrt(sum(xtra(i)^2 + ytra(i)^2 + ztra(i)^2));
    rotmat = pointreg([b_aut e1aut e2aut], [b_man e1man e2man]);
    if abs(rotmat(1,3)) == 1
        error('Angle for y-rot is pi.');
    end
    xrot(i) = atan2d(-rotmat(2,3), rotmat(3,3));
    yrot(i) = asind(rotmat(1,3));
    zrot(i) = atan2d(-rotmat(1,2), rotmat(1,1));
    magrot(i) = acosd((trace(rotmat)-1)/2);
end

%% Print stats.

fprintf('             Min    Max   Mean    Std\n\n');
stat = @(p)[min(abs(p)) max(abs(p)) mean(abs(p)) std(abs(p))];
fprintf('T(x) abs: '); fprintf('%6.2f ',stat(xtra)); fprintf('\n');
fprintf('T(y) abs: '); fprintf('%6.2f ',stat(ytra)); fprintf('\n');
fprintf('T(z) abs: '); fprintf('%6.2f ',stat(ztra)); fprintf('\n');
fprintf('R(x) abs: '); fprintf('%6.2f ',stat(xrot)); fprintf('\n');
fprintf('R(y) abs: '); fprintf('%6.2f ',stat(yrot)); fprintf('\n');
fprintf('R(z) abs: '); fprintf('%6.2f ',stat(zrot)); fprintf('\n\n');
fprintf('Magn tra: '); fprintf('%6.2f ',stat(magtra)); fprintf('\n');
fprintf('Magn rot: '); fprintf('%6.2f ',stat(magrot)); fprintf('\n\n');
fprintf('Magn eye: '); fprintf('%6.2f ',stat([mag_e1;mag_e2])); fprintf('\n\n');

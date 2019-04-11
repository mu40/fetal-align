% Remove protruding voxels (once) from each slice of a 3D mask if they are not
% included on the closest slice on at least one side.
function out = cutprot(volmask)

numvox = squeeze(sum(sum(volmask)))';
ind = find(numvox);
out = volmask;
for i = ind
    if i == ind(1)
        out(:,:,i) = volmask(:,:,i) & volmask(:,:,i+1);
    elseif i == ind(end)
        out(:,:,i) = volmask(:,:,i) & volmask(:,:,i-1);
    else
        out(:,:,i) = volmask(:,:,i) & (volmask(:,:,i-1) | volmask(:,:,i+1));
    end
end

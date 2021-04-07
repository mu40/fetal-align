% Detect MSERs on each slice of a 3D image.
function [bw,numbw,snum] = slcmser(vol, delta, minvox, maxvox)

snum = nan(10000, 1);
bw = cell(10000, 1);
numbw = 0;
sz = size(vol);
minvox = max([minvox 1]);
maxvox = min([maxvox prod(sz(1:2))]);
for i = 1:sz(3)
    im = vol(:,:,i);
    im = uint8(floor(255*im));
    im = 255 - im;
    [msers, nummser] = mser(im, delta, minvox, maxvox);
    ind = numbw + (1:nummser);
    snum(ind) = i * ones([1 nummser]);
    bw(ind) = msers;
    numbw = numbw + nummser;
end
snum = snum(1:numbw);
bw = bw(1:numbw);

% Iteratively add voxels to each slice of a 3D mask that are included on
% the closest non-zero slices on both sides.
function volmask = closegaps(volmask)

tmp = volmask;
slcind = 1:size(volmask, 3);
change = 1;
while change
    for i = slcind
        masksize = squeeze(sum(sum(tmp)))';
        j = find(slcind<i & masksize>0, 1, 'last');
        k = find(slcind>i & masksize>0, 1, 'first');
        if not(isempty(j)) && not(isempty(k))
            tmp(:,:,i) = tmp(:,:,i) | ( tmp(:,:,j) & tmp(:,:,k) );
        end
    end
    change = not(isequal(tmp, volmask));
    volmask = tmp;
end

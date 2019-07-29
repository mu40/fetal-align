function im = stretchcon(im, pctile)

if nargin() < 2
	pctile = [0.1 99.9];
end

bounds = percentile(im(~isnan(im)), pctile);
im = im-bounds(1);
if diff(bounds) > 0
    im = im / (bounds(2)-bounds(1));
end
im(im<0) = 0;
im(im>1) = 1;

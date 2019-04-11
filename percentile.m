function out = percentile(dat, pctile)

dat = squeeze(dat);

if length(pctile) ~= length(pctile(:))
    error('Percentile must be a scalar or a vector');
end
if length(dat)==length(dat(:)) && size(dat,1)==1
    dat = dat';
end

x = [0, ((0.5 : (length(dat)-0.5)) ./ length(dat)) * 100, 100];
y = [min(dat); sort(dat); max(dat)];
out = interp1(x, y, pctile);

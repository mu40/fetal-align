% Rectangular grid in N-D space. Similar to built-in ndgrid function, but
% offering more flexibility in specifying the dimension.
function out = ndarray(start, stop, step)

if nargin() < 3 || isempty(step)
    step = ones(size(start));
end

dim = fix((stop-start)./step) + 1;
numdim = numel(start);
numpoint = prod(dim);

if numel(step) == 1
    step = repmat(step, [1 numdim]);
end

if numdim == 1
    out = reshape(start:step:stop, [], 1);
    return
end

out = nan([numpoint numdim]);
for i = 1:numdim
    x = start(i):step(i):stop(i);
    x = x(1:dim(i)); % Avoid rounding problems.
    sz = ones([1 numdim]);
    sz(i) = numel(x);
    x = reshape(x, sz); % Reshape to singleton-dim array except current.
    sz = dim;
    sz(i) = 1;
    x = repmat(x, sz); % Repeat input along all other dimensions.
    out(:,i) = x(:);
end

% Assign field name/value pairs to structure if unset.
function structure = setdefault(structure, varargin)
    assert(mod(numel(varargin), 2) == 0)
    numfield = numel(varargin) / 2;
    for i = 1:numfield
        field = varargin{2*i-1};
        value = varargin{2*i};
        assert(ischar(field) && isscalar(value) && isnumeric(value))
        if isfield(structure, field)
            continue
        end
        structure.(field) = value;
    end
end

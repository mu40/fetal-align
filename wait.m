% Prompt user for key press, optionally printing a message that is
% automatically divided into multiple lines.
function wait(varargin)

linelen = 80;

msg = cellfun(@num2str, varargin, 'uniformoutput', 0);
msg = sprintf('%s ', msg{:});
if length(msg) > 1
    msg(end) = newline(); % Replace space.
    increment = linelen + 1;
    space = char(32);
    old = 1;
    pos = increment;
    while pos < length(msg)
        while msg(pos) ~= space && pos > old
            pos = pos - 1;
        end
        if pos == old
            pos = pos + increment;
            while msg(pos) ~= space && pos < length(msg)
                pos = pos + 1;
            end
        end
        if msg(pos) == space && pos < length(msg)
            msg(pos) = newline();
        end
        old = pos;
        pos = pos + increment;
    end
    fprintf(msg);
end

fprintf('Press any key to continue ...');
pause();
fprintf(' resuming\n\n');

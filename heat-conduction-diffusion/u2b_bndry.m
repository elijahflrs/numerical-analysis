function [u] = u2b_bndry(t)
    u = [t*0;                   % for t < 2
        (t*0 + 3) .* (t >= 2)]; % for t >= 2
end
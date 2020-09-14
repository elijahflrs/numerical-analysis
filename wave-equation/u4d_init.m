function [u] = u4d_init(x)
    u = (x >= 0 & x <= 1) .* x;
    u = u + (x > 1 & x <= 4) .* ((-1/3)*x+(4/3));
end

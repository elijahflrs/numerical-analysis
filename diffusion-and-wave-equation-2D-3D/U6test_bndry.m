function [U] = U6test_bndry( t, nx, ny )
    U = -Inf*ones( nx, ny );

    U([1, end], :) = 0;
    U(:, [1, end]) = 0;
end

function [U] = U6c_bndry( t, n_x, n_y )
    U = -Inf*ones( n_x, n_y );
 
    for j = 1:n_x
        for k = 1:n_y
            % (x, y) is a point in [0, 1] x [0, 1]
            x = 2*(j - 1)/(n_x - 1);
            y = (k - 1)/(n_y - 1);
        
            % Determine if a point is outside an ellipse
            if (x-1)^2/0.49 + (y- 0.5)^2/0.249 >= 1
                U(j, k) = NaN;
            end
        end
    end
end

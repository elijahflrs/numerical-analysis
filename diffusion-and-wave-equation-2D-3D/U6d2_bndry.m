function [U] = U6d2_bndry( t, nx, ny, n_z)
    U = -Inf*ones( 21, 21, 21 );
 
    for ix = 1:21
        for iy = 1:21
            for iz = 1:21
                x = (ix - 1)/(nx - 1);
                y = (iy - 1)/(ny - 1);
                z = (iz - 1)/(n_z - 1);
                dist = sqrt( (x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2 );

                if dist >= 1.1
                    U(ix, iy, iz) = NaN;
                elseif dist >= 0.5
                    U(ix, iy, iz) = sin( 1.72*t );
                end
            end
        end
    end
end

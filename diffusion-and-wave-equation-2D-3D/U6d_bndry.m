function [U] = U6d_bndry( t, n_x, n_y, n_z )
    U = -Inf*ones( n_x, n_y, n_z );
    U(end,:,:) = sin(t).*(t <= 2*pi);   % sinusoidal signal

    U(:,[1,end],:) = NaN;
    U(:,:,[1,end]) = NaN;
    
    
    for i = 1:n_x
        for j = 1:n_y
            for k = 1:n_z
                x = (i - 1)/(n_x - 1);
                y = (j - 1)/(n_y - 1);
                z = (k - 1)/(n_z - 1);


                %if the current coordinate follows a parabolic disc shape
                % with a focus of 0.25 from the vertex
                if ((z-0.5)^2 + (y-0.5)^2 >= x)
                    U(i, j, k) = NaN;
                end
            end
        end
    end

    
end

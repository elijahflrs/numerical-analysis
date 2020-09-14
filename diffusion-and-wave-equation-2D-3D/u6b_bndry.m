function [U] = u6b_bndry( t, n1, n2, n3 )
    c_bath = 1;
    T_bath = 19;
 
    if t <= T_bath
        val = c_bath;
    else
        val = NaN;
    end
    
    U = -Inf*ones( n1, n2, n3 );
    U(:, :, [1, end]) = val;
 
    for i = 1:n1
        for j = 1:n2
            for k = 1:n3
                x = (i - 1)/(n1 - 1);
                y = (j - 1)/(n2 - 1);
 
                r = sqrt( (x - 0.5)^2 + (y - 0.5)^2 );
 
                if r >= 0.5
                    U(i, j, k) = val;
                end
            end
        end
    end
end




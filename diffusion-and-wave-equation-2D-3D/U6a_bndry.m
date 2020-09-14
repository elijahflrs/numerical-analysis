function [U] = U6a_bndry( t, nx, ny )

n = nx;
U = -Inf*ones( n, n );
 
for j = 1:n
    for k = 1:n
        x = (j - 1)/(n - 1);
        y = (k - 1)/(n - 1);

        r = sqrt( (x - 0.5)^2 + (y - 0.5)^2 );

        if r < 0.1
            U(j, k) = (t <= 0.1) * 200    +    (t > 0.1 & t <= 1) * ((-1300/9) * t + (1930/9));  % inner wire at 200 for 0.1 seconds, then cools down to 70
        elseif r >= (n - 1)/(2*n)
            if y <= 0.5
                U(j, k) = 0;   % 0 on upper boundary
            else
                U(j, k) = NaN; % insulated on lower boundary
            end
        end
    end
end

end

ccc
n = 51;
U = -Inf*ones( n, n );
 
for j = 1:n
    for k = 1:n
        x = (j - 1)/(n - 1);
        y = (k - 1)/(n - 1);
        
        r = sqrt( (x - 0.5)^2 + (y - 0.5)^2 );
        
        if r < 0.1
            U(j, k) = 70;
        elseif r >= (n - 1)/(2*n)
            if y <= 0.5
                U(j, k) = 0;
            else
                U(j, k) = NaN;
            end
        end
    end
end
 
U_steady = laplace2d( U );
t_int = [0 1];
nt = 4502;
h = 0.02;

kappa = 0.43;
[t, U6a] = diffusion2d( kappa, h, U_steady, @U6a_bndry, t_int, nt );
max(U6a( round( (n + 1)/2 ), end-1, :))

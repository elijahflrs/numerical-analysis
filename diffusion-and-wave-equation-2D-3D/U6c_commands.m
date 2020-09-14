n_x = 201;
n_y = 101;
U6c_init = zeros( n_x, n_y );
dU6c_init = zeros( n_x, n_y );

for j = 1:n_x
    for k = 1:n_y
        x = 2*(j - 1)/(n_x - 1);
        y = (k - 1)/(n_y - 1);
        
        U6c_init(j, k) = exp( -1000*((x - 0.5)^2 + (y - 0.5)^2) );
    end
end

mesh(U6c_init)



[t6c, U6c] = wave2d( 1, 1, U6c_init, dU6c_init, @U6c_bndry, [0, 350], 710 );
[z_min, t_min] = min( U6c(151, 51, :) )
%t( t_min )

%frames6c = animate( U6c );

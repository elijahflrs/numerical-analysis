ccc
n_x = 33;
n_y = 33;
n_z = 33;
U6d_init = zeros( n_x, n_y, n_z );
dU6d_init = zeros( n_x, n_y, n_z );
[t6d, U6d] = wave3d( 1, 1, U6d_init, dU6d_init, @U6d_bndry, [0, 60], 150 );
%frames6d = animate( U6d, [0, 60] );



%disp(U6d(:,:,:,3))
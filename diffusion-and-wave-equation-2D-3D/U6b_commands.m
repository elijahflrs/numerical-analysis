format long
n = 21;
c_bath = 1;
U_init = zeros( n, n, n );
U_init(:, :, [1, end]) = c_bath;
 
for i = 1:n
    for j = 1:n
        for k = 1:n
            x = (i - 1)/(n - 1);
            y = (j - 1)/(n - 1);
        
            r = sqrt( (x - 0.5)^2 + (y - 0.5)^2 );
        
            if r >= 0.5
                U_init(i, j, k) = 1;
            end
        end
    end
end

% Solve for the diffusion
T_final = 19;
n_t = 500;
[t, U_out] = diffusion3d( 4, 1, U_init, @u6b_bndry, [0, T_final], n_t );

% Find the minimum and maximum concentrations in the cylinder for each point in time
c_max = zeros( 1, n_t );
c_min = zeros( 1, n_t );

for k = 1:n_t
    c_max(k) = max( max( max( U_out(:, :, :, k) ) ) );
    c_min(k) = min( min( min( U_out(:, :, :, k) ) ) );
end

% Plot the minimum and maximum temperatures
plot( t, c_max );
hold on;
plot( t, c_min );


delta = max(abs(c_max(end)-1),abs(c_min(end)-1))
%disp(10 - T_final - log2( delta ))




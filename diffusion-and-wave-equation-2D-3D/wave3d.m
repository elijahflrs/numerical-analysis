% eeflores
% wave3d.m
% This function aims to find a function u(t, x, y, z) that 
% satisfies the wave equation in three dimensions. This solution requires 
% defined initial conditions, as well as defined Dirichlet or insulated 
% conditions that will affect the form of the final solution.
%
% Parameters
% ==========
%    c - The wave constant.
%    t_int - A row vector [t_int, t_final] that is the range of time that is 
%            being solved for.
%
%    U_init - Function handle giving the initial values. 
%    dU_init - Function handle giving the initial speed of a wave.
%    U_bndry - Function handle giving the initial boundary conditions for 
%              the solution for each time in range [t_initial, t_final].
%
%    h - The interval length to discretize the domain into; i.e. the region 
%        in the spatial dimension is broken up into points separated by h.
%    n_t - The interval length to discretize the domain into; i.e. the 
%          region in the spatial dimension is broken up into points separated by h.
%
% Return Values
% =============
%    t - A vector of length n_t that contains time values in range [t_inital, t_final].
%    U_out - A matrix of size size nx x ny x x nz x nt u values approximating the solution. 

function [t, U_out] = wave3d( c, h, U_init, dU_init, U_bndry, t_int, nt )

% Argument and Error Checking
% ===========================
% if U_init or dU_init is a 2D array
if ndims(U_init) ~= 3 || ndims(dU_init) ~= 3
    throw(MException('MATLAB:invalid_argument', ...
        'The argument U_init or dU_init is not a 3D array.'));
end
% if U_bndry s a function handle
if ~isa(U_bndry, 'function_handle')
    throw(MException('MATLAB:invalid_argument', ...
        'The argument U_bndry is not a function handle.'));
end
% if t_int is a row vector of two values
if ~all(size(t_int) == [1,2])
    throw(MException('MATLAB:invalid_argument', ...
        'The argument t_int is not a row vector of 2 values.'));
end
% if n_t is a positive integer
if floor(nt) ~= nt && nt <= 0
    throw(MException('MATLAB:invalid_argument', ...
        'The argument n_t is not a positive integer.'));
end
if ~isscalar(c) || ~isscalar(h)
    throw(MException('MATLAB:invalid_argument', ...
        'the argument c or h is not a scalar value.'));
end




% Initialization
% ==============
ti = t_int(1);                  % get initial and final time values
tf = t_int(2);
dt = (tf - ti)/(nt - 1);        % get dt
t = linspace( ti, tf, nt );     % construct time vector t

[nx, ny, nz] = size( U_init ); 

U_out = zeros( nx, ny, nz, nt );    % construct the solution matrix
U_out(:, :, :, 1) = U_init;        % include the initual values 

r = (c*dt/h)^2;                 % ratio to check to ensure convergence

if r >= 0.25
    nt_suggest = ceil((c*(tf-ti)) / h + 1);
    warning( 'MATLAB:questionable_argument', ...
           'The coefficient r = (c*dt / h)^2 = %d >= 0.25. Consider using n_t = %d', ...
            r, nt_suggest)
end


% Solving
% =======
U_out(:, :, :, 2) = U_out(:, :, :, 1) + dt*dU_init;   % solve for t = t2 using Euler's method

for it = 3:nt      % for each time value
    U_out(:, :, :, it) = U_bndry( t(it), nx, ny, nz );  % find and implement the boundary condition for t(it)

    for ix = 1:nx      % for each spatial value in the x direction
        for iy = 1:ny  % for each spatial value in the y direction
            for iz = 1:nz
                if U_out(ix, iy, iz, it) == -Inf       % if the current spatial value is unknown
                    Utmp = U_out(ix, iy, iz, it - 1);  
                    U_out(ix, iy, iz, it) = 2*Utmp - U_out(ix, iy, iz, it - 2); % use the value from the time just before, solve for the current value

                    for dxyz = [-1 1 0 0 0 0; 0 0 -1 1 0 0; 0 0 0 0 1 -1]  % using the points around it (i.e. up, down, left right)
                        dix = ix + dxyz(1);
                        diy = iy + dxyz(2);
                        diz = iz + dxyz(3);


                        if ~isnan( U_out(dix, diy, diz, it - 1) ) % if it is not an insulated boundary, solve using finite difference eqn
                            U_out(ix, iy, iz, it) = U_out(ix, iy, iz, it) + ...
                                r*( U_out(dix, diy, diz, it - 1) - Utmp );
                        end
                    end
                end
            end
        end
    end
end





end

% wave1d.m
% Given an initial state, this function aims to find a function u(t,x)
% that satisfies the wave equation up to some final time t_final
% using a finite difference method. This solution requires two initial 
% boundary conditions to define the solution.
%
% Parameters
% ==========
%    c - The wave constant.
%    x_int - The interval in which the x-values (spatial values) are bounded by.
%    t_int - The interval in which your t-values (time values) are bounded by. 
%
%    u_init - Function handle giving the initial x-values in range [a, b].
%    du_init - Function handle giving the initial speed of the wave.
%    u_bndry - Function handle giving the two boundary conditions for each time in range [t0, t_final]. 
%
%    n_x - The number of points to divide the x-values (spatial values) into.
%    n_t - The number of points to divide the t-values (time values) into. 
%
% Return Values
% =============
%    x_out - A vector of x values of the n_x points in range [a,b] used in the solution.
%    t_out - A vector of t values of the n_t points in range [t_intitial, t_final] that were used in the solution.
%    U_out - A n_x x n_t matrix of u values that approximate the solution. 


function [x_out, t_out, U_out] = wave1d( c, x_int, n_x, t_int, n_t, u_init, du_init, u_bndry )

% ----- ERROR CHECKING -----
% Ensure that the input values into the function are of the correct format 
% (e.g. scalar values do not have a string as an input).
if ~isscalar( c ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a scalar' ) );
end
if ~all(size(x_int) == [1,2]) || ~all(size(t_int) == [1,2])
    throw(MException('MATLAB:invalid_argument', ...
        'the argument x_int or t_int is not a 1 x 2 row vector.'));
end
if ~isscalar(n_x) || ~isscalar(n_t)
    throw(MException('MATLAB:invalid_argument', ...
        'the argument n_x or n_t is not a scalar value.'));
end
if ~isa(u_init, 'function_handle') || ~isa(u_bndry, 'function_handle') || ~isa(du_init, 'function_handle')
    throw(MException('MATLAB:invalid_argument', ...
        'the argument u_init, du_init, or u_bndry is not a function handle.'));
end


% ----- INITIALIZATION -----
% Extract values needed to initialize the range for x and t values
a = x_int(1);
b = x_int(end);
t0 = t_int(1);
tf = t_int(end);
h = (b-a)/(n_x-1); 
dt = (tf-t0)/(n_t-1);
r = (c*dt / h)^2; % to check if the coefficient is less than 1

% if the coefficient is too large, issue a warning 
if r >= 1
    nt_suggest = ceil((c*(tf-t0)) / h + 1);
    warning( 'MATLAB:questionable_argument', ...
           'The coefficient r = (c*dt / h)^2 = %d >= 1. Consider using n_t = %d', ...
            r, nt_suggest)
end

% Construct the nx x nt matrix using zeros
U = zeros(n_x,n_t);

% Construct the initial nx values of x-values in range [a,b]. Similarly, construct the nt values of t-values in range [t0, tf].  
x_vec = linspace(a, b, n_x);
t_vec = linspace(t0, tf, n_t);


% Pass them into the function handles to determine the initial and boundary values
% Assign the initial values to the first column of the nx x nt matrix.
% Assign the boundary values in the first and last rows.U(:,1) = u_init(x_vec);
U(:,1) = u_init(x_vec)';            % initial spatial values
du_vec = du_init(x_vec)';           % initial speed of wave
boundaries = u_bndry(t_vec);
U(1,2:end) = boundaries(1,2:end);   % boundary value at x = a for all t
U(n_x,2:end) = boundaries(2,2:end); % boundary value at x = b for all t


% ----- SOLVING -----
% solve at t = t2 using Euler's method
U(:,2) =  U(:,1) + dt * du_vec;
if isnan(boundaries(1,:)) % if upper boundary is insulated
    U(1,2) = (4/3)*U(2,2) - (1/3)*U(3,2);
end
if isnan(boundaries(2,:)) % if lower boundary is insulated
    U(n_x,2) = (4/3)*U(n_x-1,2) - (1/3)*U(n_x-2,2);
end
    
% solve u values at each time step after t2
for k = 2:n_t-1 
    U(2:end-1,k+1) = 2*U(2:end-1,k) - U(2:end-1,k-1) + r*diff(U(:,k),2);  % finite difference equation
    
    if isnan(boundaries(1,:)) % if upper boundary is insulated
        U(1,k+1) = (4/3)*U(2,k+1) - (1/3)*U(3,k+1);
        
    end
    if isnan(boundaries(2,:)) % if lower boundary is insulated
        U(n_x,k+1) = (4/3)*U(n_x-1,k+1) - (1/3)*U(n_x-2,k+1);
    end
end

% output final solution
x_out = x_vec;
t_out = t_vec;
U_out = U;

end

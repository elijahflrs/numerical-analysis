% diffusion1d.m
% Given an initial state, this laboratory aims to find a function u(t,x)
% that satisfies the diffusion equation up to some final time t_final. 
% This solution requires two initial boundary conditions to define the
% solution.
%
% Parameters
% ==========
%    kappa - The diffusivity constant.
%    x_rng - The interval in which the x-values (spatial values) are bounded by.
%    t_rng - The interval in which your t-values (time values) are bounded by. 
%
%    u_init - Function handle giving the initial x-values in range [a, b].
%    u_bndry - Function handle giving the two boundary conditions for each time in range [t0, t_final]. 
%
%    nx - The number of points to divide the x-values (spatial values) into.
%    nt - The number of points to divide the t-values (time values) into. 
%
% Return Values
% =============
%    x_out - A vector of x values of the n_x points in range [a,b] used in the solution.
%    t_out - A vector of t values of the n_t points in range [t_intitial, t_final] that were used in the solution.
%    U_out - A n_x x n_t matrix of u values that approximate the solution. 


function [x_out, t_out, U_out] = diffusion1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )



% ----- ERROR CHECKING -----

% Ensure that the input values into the function are of the correct format (e.g. scalar values do not have a string as an input).
if ~isscalar( kappa ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument kappa is not a scalar' ) );
end
if ~isscalar(nx) || ~isscalar(nt)
    throw(MException('MATLAB:invalid_argument', ...
        'the argument nx or nt is not a scalar value.'));
end
if ~all(size(x_rng) == [1,2]) || ~all(size(t_rng) == [1,2])
    throw(MException('MATLAB:invalid_argument', ...
        'the argument nx or ny is not an integer.'));
end
if ~isa(u_init, 'function_handle') || ~isa(u_bndry, 'function_handle')
    throw(MException('MATLAB:invalid_argument', ...
        'the argument u_init or u_bndry is not a function handle.'));
end

% Extract values needed to initialize the range for x and t values
a = x_rng(1);
b = x_rng(2);
t0 = t_rng(1);
t_final = t_rng(2);


% Check if the coefficient is less than 0.5
h = (b-a)/(nx-1); 
coeff = kappa*((t_final-t0)/(nt-1))/ h^2;

% If it is too large, provide the smallest integer value of n_t that can be used to bring the ratio under 0.5. 
if coeff >= 0.5
    nt_suggest = ceil((kappa / (0.5*h^2)) * (t_final-t0) + 1);
    error_message = ['the ratio kappa*dt/h^2 = ', num2str(coeff), ' >= 0.5, consider using nt = ', num2str(nt_suggest)];
    throw(MException('MATLAB:invalid_argument', ...
        error_message));
end




% ----- INITIALIZATION -----

% Construct the nx x nt matrix using zeros
U = zeros(nx,nt);

% Construct the initial nx values of x-values in range [a,b]. Similarly, construct the nt values of t-values in range [t0, tfinal].  
x_vec = linspace(a, b, nx);
t_vec = linspace(t0, t_final, nt);

% Pass them into the function handles to determine the initial and boundary values. 
% Assign the initial values to the first column of the nx x nt matrix.
% Assign the boundary values in the first and last rows.
U(:,1) = u_init(x_vec);
boundaries = u_bndry(t_vec);
U(1,2:end) = boundaries(1,2:end);
U(nx,2:end) = boundaries(2,2:end); 




% ----- SOLVING -----

% Using the finite difference formula, calculate u2,2 through u_(n_x-1,2), then u2,3 through u_(n_x-1,3) and so on.
% Calculate the next rows using the formula as well.
for k = 1:nt-1
    for i = 2:nx-1
        U(i,k+1) = U(i,k) + coeff * (U(i-1,k) - 2*U(i,k) + U(i+1,k));
    end
end


% Set the final values as the output
U_out = U;
x_out = x_vec;
t_out = t_vec;


end


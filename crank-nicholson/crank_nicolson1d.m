% crank_nicolson1d.m
% Given an initial state, this function aims to find a function u(t,x)
% that satisfies the diffusion equation up to some final time t_final
% using a finite difference method This solution requires two initial 
% boundary conditions to define the solution.
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


function [x_out, t_out, U_out] = crank_nicolson1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )

% ----- ERROR CHECKING -----
% Ensure that the input values into the function are of the correct format 
% (e.g. scalar values do not have a string as an input).
if ~isscalar( kappa ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument kappa is not a scalar' ) );
end
if ~isscalar(nx) || ~isscalar(nt)
    throw(MException('MATLAB:invalid_argument', ...
        'the argument nx or nt is not a scalar value.'));
end
if  nx < 4
    throw(MException('MATLAB:invalid_argument', ...
        'the argument nx must be greater than 4.'));
end
if ~all(size(x_rng) == [1,2]) || ~all(size(t_rng) == [1,2])
    throw(MException('MATLAB:invalid_argument', ...
        'the argument nx or ny is not an integer.'));
end
if ~isa(u_init, 'function_handle') || ~isa(u_bndry, 'function_handle')
    throw(MException('MATLAB:invalid_argument', ...
        'the argument u_init or u_bndry is not a function handle.'));
end


% ----- INITIALIZATION -----
% Extract values needed to initialize the range for x and t values
a = x_rng(1);
b = x_rng(2);
t0 = t_rng(1);
t_final = t_rng(2);


% Check if the coefficient is less than 0.5
h = (b-a)/(nx-1); 
r = kappa*((t_final-t0)/(nt-1))/ h^2;

% if the coefficient is too large, issue a warning 
if r >= 0.5
    warning( 'MATLAB:questionable_argument', ...
           'the arguments of %d and %d are sub-optimal', ...
            a, b )
end

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
for k = 1:nt-1 % solve u values at each time step
    
    % ---- initialize linear equation
    known_u = U(:,k);  % isolate known u values at column k
    M = diag( (2*(1+r)) * ones(nx-2,1) ) + diag(-r * ones(nx-3,1), -1) + diag(-r * ones(nx-3,1), 1); % construct M matrix that contains coefficients
    
    % construct b matrix
    b_diff = diff(known_u,2);
    b = 2*known_u(2:end-1) + r*b_diff;
    
    % --- check if boundaries are isolated and modify matrix accordingly
    if isnan(boundaries(1,:)) & isnan(boundaries(2,:))  % if both boundaries are insulated
        M(1,1) = 2 + 2/3 * r;
        M(1,2) = -(2/3)*r;
        M(end,end) = 2 + 2/3 * r;
        M(end,end-1) = -(2/3)*r;
    elseif isnan(boundaries(1,:))                 % if left boundary is insulated
        M(1,1) = 2 + 2/3 * r;
        M(1,2) = -(2/3)*r;
        b(end) = b(end) + r*boundaries(2,k);      % add only left boundary value
    elseif isnan(boundaries(2,:))                 % if right boundary is insulated  
        M(end,end) = 2 + 2/3 * r;
        M(end,end-1) = -(2/3)*r;
        b(1) = b(1) + r*boundaries(1,k);          % add only right boundary value
    else
        b(1) = b(1) + r*boundaries(1,k);          % add both boundary values
        b(end) = b(end) + r*boundaries(2,k);
    end
    
    
    % ---- solve linear equation
    u_kplus1 = M\b;

    % ---- assign to next column of solution matrix
    U(2:end-1, k+1) = u_kplus1;
    
    if isnan(boundaries(1,:)) & isnan(boundaries(2,:))    % if both boundaries are insulated
        U(1,k+1) = (4/3)*U(2,k+1) - (1/3)*U(3,k+1);        % retreive left boundary value
        U(nx,k+1) = (4/3)*U(nx-1,k+1) - (1/3)*U(nx-2,k+1); % retreive right boundary value
    elseif isnan(boundaries(1,:))                          % if left boundary is insulated
        U(1,k+1) = (4/3)*U(2,k+1) - (1/3)*U(3,k+1);        % retreive left boundary value
    elseif isnan(boundaries(2,:))                          % if right boundary is insulated
        U(nx,k+1) = (4/3)*U(nx-1,k+1) - (1/3)*U(nx-2,k+1); % retreive right boundary value
    end
    
    
    
end

% output final solution
x_out = x_vec;
t_out = t_vec;
U_out = U;
end


%     b_diff = zeros(nx-2,1);
%     for j = 1:nx-2
%         b_diff(j,1) = known_u(j) - 2 * known_u(j+1) + known_u(j+2);
%     end

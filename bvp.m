% bvp 
%
% This function is meant to approximate a function u(x) and its derivatives
% from an LODE that equals a forcing function g(x), all while considering
% the BVP values that constrain the solution.
%
%
% Parameters
% ==========
%    c - 1x3 vector [c1, c2, c3] of constant coefficients for the LODE.
%    x_int - x-values [a,b] that bound the BVP.
%    u_int = [ua,ub], the rresult of u(a) and u(b) that bound the BVP.
%
%    g - The forcing function that the LODE is equal to.
%
%    n - The number of intervales to divide the solution into.
%
% Return Values
% =============
%    x - The n x-values used in the solution.
%    u - The n u-value solutions that were found. 

function [x, u] = bvp( c, x_int, u_int, g, n )

% Argument checking
if ~all( size( c ) == [1, 3] ) 
    throw( MException( 'MATLAB:invalid_argument', ...
    'the argument c is not a 3-dimensional row vector' ) );
end
if ~all( size( x_int ) == [1, 2] )
    throw( MException( 'MATLAB:invalid_argument', ...
    'the argument x is not a 2-dimensional row vector' ) );
end
if ~all( size( u_int ) == [1, 2] )
    throw( MException( 'MATLAB:invalid_argument', ...
    'the argument u is not a 2-dimensional row vector' ) );
end
if ~isscalar( n ) 
    throw( MException( 'MATLAB:invalid_argument', ...
    'the argument n is not a scalar' ) );
end
if ~isa( g, 'function_handle' )
    throw( MException( 'MATLAB:invalid_argument', ...
    'the argument g is not a function handle' ) );
end



% Step 1: Initialize the values u(a) = ua and u(b) = u¬b (the 
% boundaries values that the solution will be constrained on), g(x) (the 
% forcing function), [c1, c2, c3] (the constants), and n (the number of 
% intervals to divide the solution into).
a = x_int(1);
b = x_int(2);
ua = u_int(1);
ub = u_int(2);
c1 = c(1);
c2 = c(2);
c3 = c(3);

% Step 2: Create an n-sized column matrix of equally spaced x values from 
% x = a to x = b, called x. Additionally, set h = (b-a)/(n-1) to evaluate
% the solution at.
x = linspace(a, b, n)';
h = (b-a)/(n-1);

% Step 3: Using h calculate d- = (2c1 - hc2), d = (2h2c3 - 4c1), and 
% d+ = (2c1 + hc2). 
d_minus = 2*c1-h*c2;
d_ = 2*h^2*c3-4*c1;
d_plus = 2*c1+h*c2;

% Step 4: Create a (n – 2) x (n – 2) matrix with d as the diagonal 
% d- as the sub-diagonal, and d+ as the super-diagonal. Call this matrix D. 
D = diag(d_minus*ones(n-3,1), -1);
D = D + diag(d_*ones(n-2,1));
D = D + diag(d_plus*ones(n-3,1), 1);

% Step 5: Create a (n – 2) column matrix with entries of 2h2g(xk) with 
% k = 2 to k = n – 2 (k referring to the kth x-value in matrix x). To the 
% first entry, subtract d-ua. To the last entry, subtract d+ub. Call this 
% matrix G.

G = 2*h^2*g(x(2:end-1));
G(1) = G(1) - d_minus*ua;
G(end) = G(end) - d_plus*ub;

% Step 6: Solve U = D\G. 
U = D\G;

% Step 7: Create final matrix u, which are u-values used (including 
% the initial boundary conditions). The outputted variables will be matrix 
% x and u.
u = [ua; U; ub];


end
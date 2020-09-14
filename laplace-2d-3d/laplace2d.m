% laplace2d
% Given initial boundary conditions, this function aims to approximate a 
% function u in which its second derivatives in all directions (i.e. the x- 
% and y- directions) are equal to zero. 
%
% Parameters
% ==========
%    U - A n_x by n_y empty matrix that will hold the value of all points 
%        of the Laplace function that is being approximated.
%
% Return Values
% =============
%    U_out - A n_x by n_y matrix that holds the calculated solutions of all 
%            points of the Laplace function that is being approximated.

function [U_out] = laplace2d( U )

% Error and Warning Checking
% ==========================
% check if U matrix is 2D
if ndims(U) ~= 2
    throw(MException('MATLAB:invalid_argument', ...
        'The argument U is not a 2D matrix.'));
end

% check if none of the points on the edge of the array are -inf
if sum(U(1,:) == -inf) > 0 || sum(U(end,:) == -inf) > 0 || sum(U(:,1) == -inf) > 0 || sum(U(:,end) == -inf) > 0
    throw(MException('MATLAB:invalid_argument', ...
        'The argument U does not have a specified boundary (fixed or insulated) around the entire array (i.e. there is an -Inf value on the edge of the array).'));
end



% Initialization
% ==============
[n_x, n_y] = size(U);
U_out = U;               % initialize matrix U_out using U



% Mapping the unknown points to a unique number from 1 to m
% =========================================================
u_to_w = zeros( n_x, n_y );
w_to_u = zeros( 2, n_x * n_y );
m = 0;
for ix = 1:n_x
    for iy = 1:n_y
        if U(ix, iy) == -Inf
            m = m + 1;
            u_to_w(ix, iy) = m;             % for each unknown, map it to a unique number from 1 to m
            w_to_u(:, m) = [ix, iy]';       % for each w-space value, indicate the u-coordinate that it belongs to
        end
    end
end



% Creating and solving a system of linear equations
% =================================================
M = spalloc(m, m, 5*m);
b = zeros(m, 1);

for k = 1:m
    c = w_to_u(:,k);  % get u-coordinates of the unknown point we are investigating
    coordinates = [[-1 0]', [1 0]', [0 -1]', [0 1]'];   % list of the relative coordinates of the upper, lower, left, and right neighbouring point (respectively)
    for coord = 1:4
        p = c + [coordinates(1,coord), coordinates(2,coord)]';  % get u-coordinate of the neighbouring point being investigated
        i = k;                      % change to variable i for consistency with instructions
        if U(p(1), p(2)) == -Inf    % if neighbouring point is an unknown                
            j = u_to_w(p(1), p(2)); % find w-space coordinate of the neighbouring point
            M(i,i) = M(i,i) - 1;    % subtract 1 off ith diagonal of M
            M(i,j) = M(i,j) + 1;    % add 1 to (i,j)th point on M  
        elseif isnumeric(U(p(1), p(2))) && ~isnan(U(p(1), p(2)))  % if neighbouring point is a Dirichlet boundary
            M(i,i) = M(i,i) - 1;            % subtract 1 off ith diagonal of M
            b(i) = b(i) - U(p(1), p(2));    % subtract the value from ith entry of b matrix
        end
    end 
end

w = M \ b;  % solve for the unknown points



% Substituting the values back into the matrix U_out
% ===================================================
for k = 1:m
    c = w_to_u(:,k);           % get u-coordinate of the solved w value
    U_out(c(1), c(2)) = w(k);  % copy value from w into U_out
end


end
% laplace3d
% Given initial boundary conditions, this function aims to approximate a 
% function u in which its second derivatives in all directions (i.e. the x-, 
% y-, and z- directions) are equal to zero. 
%
% Parameters
% ==========
%    U - A n_x by n_y by n_z empty matrix that will hold the value of all points 
%        of the Laplace function that is being approximated.
%
% Return Values
% =============
%    U_out - A n_x by n_y by n_z matrix that holds the calculated solutions of all 
%            points of the Laplace function that is being approximated.

function [U_out] = laplace3d( U )

% Error and Warning Checking
% ==========================
% check if U matrix is 2D
if ndims(U) ~= 3
    throw(MException('MATLAB:invalid_argument', ...
        'The argument U is not a 3D matrix.'));
end

% check if none of the points on the edge of the array are -inf
if sum(sum(U(1,:,:) == -inf)) > 0 || sum(sum(U(end,:,:) == -inf)) > 0 || sum(sum(U(:,:,1) == -inf)) > 0 || sum(sum(U(:,:,end) == -inf)) > 0 || sum(sum(U(:,end,:) == -inf)) > 0 || sum(sum(U(:,1,:) == -inf)) > 0
    throw(MException('MATLAB:invalid_argument', ...
        'The argument U does not have a specified boundary (fixed or insulated) around the entire matrix (i.e. there is an -Inf value on the edge of the matrix).'));
end



% Initialization
% ==============
[n_x, n_y, n_z] = size(U);
U_out = U;               % initialize matrix U_out using U



% Mapping the unknown points to a unique number from 1 to m
% =========================================================
u_to_w = zeros( n_x, n_y , n_z);
w_to_u = zeros( 3, n_x * n_y * n_z);
m = 0;
for ix = 1:n_x
    for iy = 1:n_y
        for iz = 1:n_z
            if U(ix, iy, iz) == -Inf
                m = m + 1;
                u_to_w(ix, iy, iz) = m;             % for each unknown, map it to a unique number from 1 to m
                w_to_u(:, m) = [ix, iy, iz]';       % for each w-space value, indicate the u-coordinate that it belongs to
            end
        end
    end
end



% Creating and solving a system of linear equations
% =================================================
M = spalloc(m, m, 5*m);
b = zeros(m, 1);

for k = 1:m
    c = w_to_u(:,k);  % get u-coordinates of the unknown point we are investigating
    coordinates = [[-1 0 0]', [1 0 0]', [0 -1 0]', [0 1 0]', [0 0 -1]', [0 0 1]'];   % list of the relative coordinates of the upper, lower, left, right, front, and back neighbouring points (respectively)
    for coord = 1:6
        p = c + [coordinates(1,coord), coordinates(2,coord), coordinates(3,coord)]';  % get u-coordinate of the neighbouring point being investigated
        i = k;                            % change to variable i for consistency with instructions
        if U(p(1), p(2), p(3)) == -Inf    % if neighbouring point is an unknown                
            j = u_to_w(p(1), p(2), p(3)); % find w-space coordinate of the neighbouring point
            M(i,i) = M(i,i) - 1;          % subtract 1 off ith diagonal of M
            M(i,j) = M(i,j) + 1;          % add 1 to (i,j)th point on M  
        elseif isnumeric(U(p(1), p(2), p(3))) && ~isnan(U(p(1), p(2), p(3)))  % if neighbouring point is a Dirichlet boundary
            M(i,i) = M(i,i) - 1;          % subtract 1 off ith diagonal of M
            b(i) = b(i) - U(p(1), p(2), p(3));    % subtract the value from ith entry of b matrix
        end
    end 
end

w = M \ b;  % solve for the unknown points



% Substituting the values back into the matrix U_out
% ===================================================
for k = 1:m
    c = w_to_u(:,k);           % get u-coordinate of the solved w value
    U_out(c(1), c(2), c(3)) = w(k);  % copy value from w into U_out
end


end
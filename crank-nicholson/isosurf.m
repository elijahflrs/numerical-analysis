% isosurf( U, n )
%
% Plot isosurfaces for solutions to various equations
% in three dimensions.
%
%  Author:  Douglas Wilhelm Harder
%  Copyright (c) 2010-11 by Douglas Wilhelm Harder.  All rights reserved.
%
% The input is an m1 x m2 x m3 array A which holds
% values within the region and plots n + 1 isosurfaces
% interpolated from the smallest to the largest value.
%
% The colour spectrum is one of:
%      'rgb'     red-green-blue
%      'rb'      red-blue
%      'gb'      green-blue
%
% The transparancy can be adjusted by using the alpha( x ) command
% where the argument 0 <= x <= 1 goes from transparent to opaque.

function [] = isosurf( U, n, rng, spectrum )
    if nargin == 1
        throw( MException( 'MATLAB:invalid_argument', ...
               'expecting 2 or 3 arguments but got %d', nargin ) );
    end
    
    if nargin == 2 || (nargin >= 3 && all( size( rng ) == [0, 0] ))
    	lo = min( min( min( U ) ) );
    	hi = max( max( max( U ) ) );

    	vals = linspace( lo, hi, n );
    elseif nargin >= 3
        if ~all( size( rng ) == [1, 2] )
            throw( MException( 'MATLAB:invalid_argument', ...
                   'the range must be a 2-D row vector [a, b]' ) );
        end
        
        lo = rng(1);
        hi = rng(2);
        
        if lo >= hi
            throw( MException( 'MATLAB:invalid_argument', ...
                   'the range must be ordered but %f >= %f', lo, hi ) );
        end
        
        vals = linspace( lo, hi, n );
 
        lo = min( min( min( U ) ) );
    	hi = max( max( max( U ) ) );
    end
    
    if nargin <= 3
        spectrum = 'rgb';
    end
    
    newplot;

    for i = 1:n
        if vals(i) >= lo && vals(i) <= hi
            r = (i - 1)/(n - 1);
            
            if strcmp( spectrum, 'rgb' )
                colour = rgb( r );
            elseif strcmp( spectrum, 'rb' )
                colour = rb( r );
            elseif strcmp( spectrum, 'gb' )
                colour = gb( r );
            elseif isa( spectrum, 'function_handle' )
                colour = spectrum( r );
                
                if all( size( colour ) ~= [1, 3] )
                    throw( MException( 'MATLAB:invalid_colour', ...
                       'the colour should be a triplet [r, g, b]' ) );
                end
            else
                throw( MException( 'MATLAB:invalid_argument', ...
                   'the spectrum should be one of ''quadratic'', ''linear'', ''rb'', ''gb'' or a function handle' ) );
            end
            
            patch(                        ...
    			isosurface( U, vals(i) ), ...
	    		'FaceColor', colour,      ...
    			'EdgeColor', 'none'       ...
	    	);
        end
    end

    view( 3 );
	alpha(1/n);
	axis equal;
end

function v = rgb( r )
    v = [                                                 ...
        (-3 + 6*r)*(r >= 1/2 & r <= 2/3) + (r > 2/3),     ...
        (   10*r/3)*(r <= 3/10) + (r > 3/10 & r <= 2/3) + ...
        (3 - 3*r)*(r > 2/3),                              ...
        (r <= 3/10) + (5/2 - 5*r)*(r > 3/10 & r <= 1/2)   ...
    ];
end

function v = rb( r )
    v = [                                ...
        (2*r)*(r <= 1/2) + (r > 1/2),    ...
        0,                               ...
        (r <= 1/2) + (2 - 2*r)*(r > 1/2) ...
    ];
end

function v = gb( r )
    v = [                                ...
        0,                               ...
        (2*r)*(r <= 1/2) + (r > 1/2),    ...
        (r <= 1/2) + (2 - 2*r)*(r > 1/2) ...
    ];
end
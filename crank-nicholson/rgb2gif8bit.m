% rgb2gif8bit
%
%  Author:  Douglas Wilhelm Harder
%  Copyright (c) 2010-11 by Douglas Wilhelm Harder.  All rights reserved.
%
% Convert a 24-bit colour image into an 8-bit colour image using
% a colour map.
%
% 3^3 + 6^3 = 243
% 
% The first 3^3 = 27 points are the corners when the 256 x 256 x 256
% points are broken into 2^3 = 8 cubic regions.
%
% For the other 6^3 points, the grid of 256 x 256 x 256
% points are broken into 6^3 = 216 cubic regions where each
% region represented by the mid-point.

function [Gif8bit, color_map] = rgb2gif8bit( ImageRGB )
    color_map = zeros( 243, 3 );

    % Create the colour map
    
    % First, divide the region into 8 equally sized cubes
    % and store the corners of those cubes as the first 27
    % colours in the colour map.
    
    for r = 0:2
        for g = 0:2
            for b = 0:2
                color_map(1 + r + 3*g + 9*b, :) = [r, g, b]/2;
            end
        end
    end

    for r = 0:5
        for g = 0:5
            for b = 0:5
                color_map(28 + r + 6*g + 36*b, :) = [r, g, b]/6 + 1/12;
            end
        end
    end

    % Next, for each triplet of RGB colours in the 24-bit image,
    % replace it with the colour in the colour map that 
    
    m = size( ImageRGB, 1 );
    n = size( ImageRGB, 2 );
    Gif8bit = zeros( m, n );
    
    for i = 1:m
        for j = 1:n
            Gif8bit(i, j) = rgb2map( double( reshape( ImageRGB(i, j, :), 1, 3 ) ) );
        end
    end
end

function [idx] = rgb2map( rgb )
    % The background of Matlab plots are [204, 204, 204]
    % It looks cleanest to display this one specific colour as white.
    if all( rgb == 204 )
        idx = 27;
        return;
    end
    
    % Find the closest triplets for the values 
    map1 = round( (2/255)*rgb );               % In {0, 1, 2}^3
    map2 = min( 5, floor( (6/255)*rgb ));      % In {0, 1, 2, ..., 5}^3
    
    % Map back to the closest color in [0, 255]^3
    rgb1 = (255/2)*map1;
    rgb2 = (255/6)*(map2 + 0.5);
    
    % Choose which ever colour is closer
    if norm( rgb - rgb1 ) < norm( rgb - rgb2 )
        idx = 1 + map1(1) + map1(2)*3 + map1(3)*9;
    else
        idx = 28 + map2(1) + map2(2)*6 + map2(3)*36;
    end
end

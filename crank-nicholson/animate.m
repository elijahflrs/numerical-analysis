% anmiate
%
% animate( U ) plots an animation of the slices U(:,:,i) or U(:,:,:,i)
% of the argument U
%
% animate( U, [a, b] ) restricts the range of the plot to [a, b]
%
% animate( U, [a, b], filename.gif ) or animate( U, [], filename.gif )
% stores the movie to the GIF file with the name filename.gif.
%
% You CANNOT block the Plot window while the images are being captured.
% MATLAB will capture the image that is visible on the screen in front
% of the Plot window. 

function [frames] = animate( U, range )
    if nargin == 1
        range = [];
    elseif nargin == 2
        if ~all( size( range ) == [1, 2] ) && ~all( size( range ) == [0, 0] ) 
            throw( MException( 'MATLAB:invalid_argument',               ...
               'the 2nd argument must be either [] or a 2-D row vector' ...
	        ) );
        end
    else
        throw( MException( 'MATLAB:invalid_argument',              ...
               'expecting 1 or 2 arguments, but got %d', nargin ...
        ) );
    end
    
    % Get some information about the array U
    dims = size( U );
    n_t = dims( end );
    ndims = length( dims );

    % Create an array of n_t cdata-colormap data structures
    frames(n_t) = struct( 'cdata', [], 'colormap', [] );

    if ndims == 2
        % If no range is specified by the user, find the minimum and
        % maximum values that appear in the array U.
        
        if ~all( size( range ) == [1, 2] )
            range = [min(min(U)), max(max(U))];
        end
        
        if range(1) == range(2)
            range = range + [-1, 1];
        end

        % If each time slice is one-dimensional, simply plot the
        % points and rescale the y-axis to encompass all possible values
        % that appear in the array U unless overridden by the user.
        
        for it = 1:n_t
            % Create a plot of the it_th time slice
            newplot;
            plot( U(:,it), 'r-' );
            ylim( range );
            title( 'dwharder' );

            % Capture the current plot as the it'th frame of the movie
            frames(it) = getframe;
        end
    elseif ndims == 3
        % If no range is specified by the user, find the minimum and
        % maximum values that appear in the array U.
        
        if ~all( size( range ) == [1, 2] )
            range = [min(min(min(U))), max(max(max(U)))];
        end
        
        if range(1) == range(2)
            range = range + [-1, 1];
        end
        
        % If each time slice is two-dimensional, use 'mesh' to
        % plot the points and rescale the z-axis to capture all possible
        % values that appear in the array U unless overridden by the user.
 
        for it = 1:n_t
            newplot;
            mesh( U(:,:,it) );
            zlim( range );
            title( 'dwharder' );

            % Capture the current plot as the it'th frame of the movie
            frames(it) = getframe;
        end
    elseif ndims == 4
        % If no range is specified by the user, find the minimum and
        % maximum values that appear in the array U.
        
        if ~all( size( range ) == [1, 2] )
            range = [min(min(min(min(U)))), max(max(max(max(U))))];
        end
        
        if range(1) == range(2)
            range = range + [-1, 1];
        end
        
        for it = 1:n_t
            % Create a plot of the it_th time slice
            newplot;

            % Change the colour range:
            %     'rgb'    red-green-blue
            %     'rb'     red-blue
            %     'gb'     green-blue
            isosurf( U(:,:,:,it), 51, range, 'rgb' );
            % Reduce the transparency by increasing the alpha value
            % alpha( 0.2 );
            
            xlim( [1, dims(1)] );
            ylim( [1, dims(2)] );
            zlim( [1, dims(3)] );
            title( 'dwharder' );
            
            % Capture the current plot as the it'th frame of the movie
            frames(it) = getframe;
        end
    else
        throw( MException( 'MATLAB:invalid_argument',                    ...
               'the first argument must be either a 2-, 3- or 4-D array' ...
        ) );
    end
end
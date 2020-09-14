% frames2gif( frames, filename, ... )
%
%  Author:  Douglas Wilhelm Harder
%  Copyright (c) 2010-11 by Douglas Wilhelm Harder.  All rights reserved.
%
% Save a sequence of frames as an animated gif.
%
% General format taken movie2gif.m by Nicolae Cindea; however, there are
% significant modifications.
%                                 http://www.iecn.u-nancy.fr/~cindea/
%
% 1. Error Checking
% Check that the 2nd argument is a string and append '.gif' if the
% file name does not already end in '.gif'.
%
% 2. Display a Progress Bar
% This will let the user see the progress of the generation of the gif.
% If any exceptions are thrown after this point, catch the exception,
% close the progress bar dialog and rethrow the excpetion.
%
% 3. Write the First Frame
% Create a new .gif file with frame converted to an image with an
% 8-bit colour map.  Set the gif to loop (an infinite loop count).
%
% 4. Append Subsequent Frames
% For each subsequent frame:
%   a. Update the progress bar,
%   b. Convert the frame to an image with an 8-bit colour map, and
%   c. Append the image to the gif file.
% 
% 5. Clean Up
% Display a full progress bar and close the dialog.
% Seeing the progress bar go to 100 % is more appropriate than having it
% at 95 % full and the closing.  Closing it early makes it appear as
% if there may have been a failure.

function [] = frames2gif( frames, filename, varargin )
    % Check the the 2nd argument is a string
    if ~ischar( filename )
        throw( MException( 'MATLAB:invalid_argument',            ...
              'the 2nd argument, ''filename'', must be a string' ...
        ) );
    end
    
    % If the last four characters of the file name are not '.gif',
    % append thouse characters.
    if length( filename ) < 4 || ~all( filename((end - 3):end) == '.gif')
        filename = [filename, '.gif'];
    end
    
    % Display a progress bar in a dialog.
    progress = waitbar( 0, sprintf( 'Saving Animation to ''%s''', filename ) ); 

    % If any exception is thrown, close the progress bar before
    % retrowing the exception
    try
        nframes = size( frames, 2 );

        % Write the first frame to a new gif file
        [RGB, ~] = frame2im( frames(1) );
        [IND, colormap] = rgb2gif8bit( RGB );
        imwrite( IND, colormap, filename, 'gif', 'LoopCount', Inf, varargin{:} );

        for k = 2:nframes
            % Update the progress bar
            waitbar( (k - 1)/nframes, progress );

            % Convert the frame into an image (a m x n x 3 array of RGB
            % colours)
            [RGB, ~] = frame2im( frames(k) );
            % A gif uses only 8 bits for colour; therefore, the 24-bit colour
            % must be converted to 2^8 = 256 colours using a colour map.
            [IND, colormap] = rgb2gif8bit( RGB );
            warning off;
            % Append the next frame to the gif file
            imwrite( IND, colormap, filename, 'gif', 'WriteMode', 'append', ...
                     'DelayTime', 0, varargin{:} );
            warning on;
        end
    catch err;
        % If there is an exception caught, close the progress bar
        % and rethrow the excetion.
        close( progress );
        rethrow( err );
    end
    
    % Complete the progress bar and the close it
    waitbar( 1, progress );
    close( progress );
end

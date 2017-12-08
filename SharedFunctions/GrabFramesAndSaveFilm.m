function [ Frame ] = GrabFramesAndSaveFilm( LpNo,Frame,isOctave )
%GrabFramesAndSaveFilm: Put inside loop to create frames for exporting a
%                   Movie. Works for Octave and MATLAB (Octave is a bit
%                   tricky so we just save images to the current dir). Note
%                   the inputs and outputs to the func are not used in
%                   Octave.
% usage #1:
% [ Frame ] = GrabFramesAndSaveFilm( LpNo,Frame,isOctave )
%
%
% Arguments: (input)
% LpNo              - The current loop number. 
%
%
% Frame             - For MATLAB's Movie making func
%
% isOctave          - Flag, if its 1 its assumed octave is being used and
%                     pictures are produced in the current dir (named
%                     monotomically).
%
% Arguments: (output)
% Frame             - For MATLAB's Movie making func
%
% Example usage:
% 
% x = 0:pi/100:2*pi;
% y = sin(x);
% plot(x,y); hold on
% Frame=getframe(gcf); %Start array 'Frame'
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
% fps=24;
% for i=1:15; 
%   x=x+pi/i
%   plot(x,y)
%   [ Frame ] = GrabFramesAndSaveFilm( i,Frame,isOctave );
% end
% if  isOctave==1
%     disp('Octave has issues Animating, use the created images and create a gif/movie online')
% elseif isOctave==0
%     % Exporting video if wanted
%     %movie(Frame,loops,fps); %movie(Frame,n(loops),fps(default12))
%     v=VideoWriter('Animation','MPEG-4');
%     v.FrameRate=fps;
%     open(v)
%     writeVideo(v,Frame)
%     disp('duration will be this many seconds, change framerate to increase')
%     disp(v.Duration)
%     close(v)
% end
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


if  isOctave==1
    
    f = warndlg('Octave is about to spit loads of images in your current directory, do you really want to continue?, ctrl+C in cmd window to quit, if you do continue make sure your in a directory that you do not mind filling with images');
    drawnow     % Necessary to print the message
        
    print 'Image'.jpg; %creating image
    oldname='Image.jpg'; %its current name 
    name=(['Image' num2str(LpNo)]); %adding number
    movefile(oldname, [name '.jpg']); %renaming
    
    
elseif isOctave==0
    
    % Store Each frame (Turn on when creating Movie)
    Frame(LpNo)=getframe(gcf); % leaving gcf out crops the frame in the movie.
    
end

end


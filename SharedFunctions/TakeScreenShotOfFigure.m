function TakeScreenShotOfFigure(num)
% TakeScreenShotOfFigure: Takes a screenshot of the currently selected 
%                   figure (gcf). Exports as a png in the current directory
%                   called 'Test(num).png'. Works well when linked with:
%                   ChangeFontSizes.m as this function inflates the fig to
%                   full screen mode before saving it.
%   
% usage #1:
% TakeScreenShotOfFigure(num)
%
% Arguments: (input)
% num              - The figure number saved when exported. 
%
% Example usage:
%
% x = 0:pi/100:2*pi;
% y = sin(x);
% plot(x,y)
% ChangeFontSizes(12,14)
% TakeScreenShotOfFigure(1)
% close
% %you will end up with a png called 'test1'
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Getting figure and inflating this to full screen size
set(gcf, 'Position', get(0,'Screensize'));
%Saving it in current dir
print(strcat('Test',num2str(num)),'-dpng')


end

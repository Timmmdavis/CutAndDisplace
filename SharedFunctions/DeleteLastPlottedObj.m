function DeleteLastPlottedObj
% DeleteLastPlottedObj: Deletes the last object plotted on a figure. For
%                   example if a line has been plotted by mistake this
%                   makes this easy to remove. Just call in cmd window. 
%               
% usage #1: 
% DeleteLastPlottedObj
%
% Arguments: (input)
% N/A
%
% Arguments: (output)
% N/A
%
% Example usage 1: Plot lines in order (RGB) then remove the last line
% % (green). 
% figure;
% plot(sind(linspace(0,360,360)),'r')
% hold on
% plot(cosd(linspace(0,360,360)),'b')
% plot(-cosd(linspace(0,360,360)),'g')
% title('3 Lines plotted, RGB')
% CloneFigure %Clone so you can see the difference
% title('3 Lines plotted, RGB, G then removed')
% DeleteLastPlottedObj
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

items = get(gca, 'Children');
delete(items(1));

end


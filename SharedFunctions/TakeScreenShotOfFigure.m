function TakeScreenShotOfFigure(num)
%Takes a screenshot of the currently selected figure (gcf)
%Exports as a png in the current directory called 'Test(num).png'

%num is the figure number you want

%   Copyright 2017, Tim Davis, The University of Aberdeen

% if you have a figure showing call like:
%TakeScreenShotOfFigure(1)
%you will end up with a png called 'test1'

%grabbing and printing at fullscreensize
set(gcf, 'Position', get(0,'Screensize')); 
print(strcat('Test',num2str(num)),'-dpng')
end

function [ FigaxisOb,FigHand ] = CloneFigure
% CloneFigure: Clones the current figure (selected) into a new plot or
%                   returns the axis handles to this and it's embedded
%                   objects for putting into a subplot. 
%               
% usage #1: Simply creates a close of the current figure
% CloneFigure 
%
% usage #2: Returns graphic objects that correspond to the current fig. 
% [ FigaxisOb,FigHand ] = CloneFigure
%
% Arguments: (input)
% N/A
%
% Arguments: (output)
% FigaxisOb         - Objects that the figure contains. 
% 
% FigHand           - Handle to the figure. 
%
% Example usage 1: Clones into new figure
% close all
% peaks; %Create fig
% pause(0.2)
% CloneFigure
%
% Example usage 2: Place current figure into subplot
% close all
% peaks; %Create fig
% pause(0.2)
% [FigaxisOb,FigHand] = CloneFigure
% Fig2=figure;
% s1=subplot(1,2,1);
% pos1=get(s1,'Position');
% delete(s1); 
% FigaxisOb=copyobj(FigaxisOb,Fig2);
% set(FigaxisOb, 'Position', pos1);
% %And if you want to delete that first fig:
% close(FigHand);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


if nargout==0 
    %We just create the figure
    copyobj(gcf,0);
else
    %Provide axis outputs
    %Axis objects
    FigaxisOb=gca;
    %The actual figure handle
    FigHand=gcf;
end


end


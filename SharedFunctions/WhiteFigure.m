function [] = WhiteFigure
% WhiteFigure: Makes the current figure turn white (both background and
%              outline).
%               
% usage #1:
% WhiteFigure
%
% Arguments: (input)
% N/A
%
% Arguments: (output)
% N/A
%
% Example usage 1:
%
% % Making MATLABS peaks plot example look better:
% % Draws a nice surface with title etc
% peaks;
% % We copy this to new figure:
% f2=copyobj(gcf,0);
% %Apply the function
% WhiteFigure
% %Now compare the two figures!
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

set(gcf,'color',[1 1 1])%gcf=Get Current Figure, to set border/figure color
set(gca,'color',[1 1 1])%gca=Get Current Axes, to set axes color
set(gca,'box','on')

end


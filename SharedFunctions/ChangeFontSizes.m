function  ChangeFontSizes( AxisFontSize,TitleFontSize )
% ChangeFontSizes: Change Font sizes of the current selected figure (either
%                  by gcf or just clicking on it). Also makes the figure
%                  background white. Effects legend entries etc also. 
%               
% usage #1:
% ChangeFontSizes( AxisFontSize,TitleFontSize )
%
% Arguments: (input)
% AxisFontSize 		- The font size of the axis labels.
%
% TitleFontSize     - The font size of the titles labels.
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
% ChangeFontSizes( 16,20 )
% %Now compare the two figures!
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Grabbing all the text in the figure. 
set(findall(gcf,'type','axes'),'fontsize',AxisFontSize)
set(findall(gcf,'type','text'),'fontSize',TitleFontSize) 
%Setting the background white
set(gcf,'color','w'); 

%changing colorbar font only if it exists
figure_handle=gcf;
cbar_handle = findobj(figure_handle,'tag','Colorbar');
if ~sum(size(cbar_handle))==0
    handle=colorbar;
    set(handle,'fontsize',AxisFontSize);
end


end


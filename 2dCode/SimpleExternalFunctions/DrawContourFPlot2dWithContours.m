function DrawContourFPlot2dWithContours( X,Y,Data,Lvlstp,cmap )
% DrawContourFPlot2dWithContours: Draws a filled contour plots with
% 					labelled contour lines shown as black. 
%                   
%               
% usage #1:
% DrawContourFPlot2dWithContours( X,Y,Data,Lvlstp,cmap )
%
% Arguments: (input)
%     X,Y          - X and Y locations of the data defined in the other
%                   input arguments. Must be on a grid.
%
%     Data 		   - The Z level or contour surface, must be on a grid and
%					have the same size as XY. 
%
% 	Lvlstp		   - The step between the contour levels that you want to
%                   have text on. Note this will look better if there are
%                   not too many decimal places.
%
%     cmap         - A colourmap that MATLAB can use. See func
%                   "colormap_cpt.m" to produce one. 
%
% Arguments: (output)
% N/A
%
%
% Example usage:
% cmap2 = colormap_cpt('Ccool-warm2');
% [x,y] = meshgrid(-2:.1:2);                                
% z = x + (-x.^2 - y.^2);
% NoOfSteps=10;
% Lvlstp=(max(z(:))-min(z(:)))/NoOfSteps; %(Alternativly just choose a
% number). 
% DrawContourFPlot2dWithContours( x,y,z,Lvlstp,cmap2 )
% figure;surf(x,y,z);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



%Creating specific contour levels, See MATLABS 
%''highlight-specific-contour-levels'' Doc. 
minD = floor(min(Data(isfinite(Data(:)))));
maxD = ceil(max(Data(isfinite(Data(:)))));
indx = minD:Lvlstp:maxD;

%Setting axis limits so supurious values do not change the overall trend. 
caxiss=[min(Data(:)),max(Data(:))];
%Rounding values, rounds the axis values to the set step.
rnd=1/Lvlstp;
caxiss = (round(caxiss*rnd))/rnd; 

%Now drawing
[~,h] =contourf(X,Y,Data,indx);
hold on
%Adding contours
contour(X,Y,Data,indx,'k');
set (h, 'levelstep', Lvlstp);
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
%Colourbars and map
colorbar;
caxis(caxiss);
if isempty(cmap)
    colormap(gray);
else
    colormap(cmap);
end 
%Labels
xlabel('x');
ylabel('y');
axis('equal'); 

%Make figure white
WhiteFigure

end


function DrawContourFPlots2d( X,Y,cmap, varargin )
% DrawContourFPlots2d: Draws multiple filled contour plots of various
%                   arguments. XY data must be on a grid as well as the
%                   arguments in varargin. If not see
%                   "DrawScatterPlots2d.m". 
%
%                   Pass in X and Y,the colourmap and various data at the
%                   XY positions and this plots each of the varargin for
%                   you. Name your input arguments well as these will be
%                   the title.
%               
% usage #1:
% DrawContourFPlots2d( X,Y,cmap,varargin )
%
% usage #2: No colourmap defined, defaults to grey. 
% DrawContourFPlots2d( X,Y,[],varargin )
%
% Arguments: (input)
%     X,Y          - X and Y locations of the data defined in the other
%                   input arguments. Must be on a grid.
%
%     cmap         - A colourmap that MATLAB can use. See func
%                   "colormap_cpt.m" to produce one. 
%
%     varargin     - Arguments you want to draw. For example the stress
%                   tensor:Sxx, Syy, Sxy. Make sure these are named well as
%                   these will be the plot title. Must be on a grid and the
%                   same size as XY.
%
%
% Example usage:
%
%  DrawContourFPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

for i=1:nargin-3
    VarName=inputname(i+3); %grabbing the input variables name, is a string
    Data=varargin{i};
    figure;
    contourf(X,Y,Data);
    if isempty(cmap)
        colormap(gray);
    else
        colormap(cmap);
    end 
    xlabel('x');
    ylabel('y');
    axis('equal'); 
    title(VarName);
    DivergingCentre( Data )
    %Make figure white
    WhiteFigure
end


function DrawScatterPlots3d( X,Y,Z,cmap, varargin )
% DrawScatterPlots3d: Draws multiple scatter plots of various
%                   arguments. XYZ data must match the size of the varargin.
%
%                   Pass in X, Y and Z the colourmap and various data at the
%                   XYZ positions and this plots each of the varargin for
%                   you. Name your input arguments well as these will be
%                   the title.
%               
% usage #1:
% DrawScatterPlots3d( X,Y,Z,cmap, varargin )
%
% usage #2: No colourmap defined, defaults to grey. 
% DrawScatterPlots3d( X,Y,Z,[], varargin )
%
% Arguments: (input)
%     X,Y,Z        - X, Y and Z locations of the data defined in the other
%                   input arguments. Must be a vector.
%
%     cmap         - A colourmap that MATLAB can use. See func
%                   "colormap_cpt.m" to produce one. 
%
%     varargin     - Arguments you want to draw. For example the stress
%                   tensor:Sxx, Syy, Sxy. Make sure these are named well as
%                   these will be the plot title. Must be a vector and
%                   match size of XY and Z.
%
%
% Example usage:
%
%  DrawScatterPlots3d( X,Y,Z,cmap2,Sxy(:),Syy(:),Sxx(:),S1(:),S2(:) );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


for i=1:nargin-4
    VarName=inputname(i+4); %grabbing the input variables name, is a string
    Data=varargin{i};
    figure;
	scatter3(X(:),Y(:),Z(:),15,Data(:));
    if isempty(cmap); colormap(gray); else; colormap(cmap); end 
    xlabel('x'); ylabel('y'); axis('equal'); 
    title(VarName);
    DivergingCentre( Data )
end


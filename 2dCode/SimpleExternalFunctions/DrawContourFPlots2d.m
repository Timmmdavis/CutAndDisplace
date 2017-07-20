function DrawContourFPlots2d( X,Y,cmap, varargin )
%DrawScatterPlots2d Pass in X Y the colourmap and data and this plots each
%of the varargin for you. 
%X and Y are the X and Y Data.
%Cmap is the colourmap you will use
% I presume your variables are named well, the title takes the input arguments
% name


%   Copyright 2017, Tim Davis, The University of Aberdeen
for i=1:nargin-3
    VarName=inputname(i+3); %grabbing the input variables name, is a string
    Data=varargin{i};
    figure;contourf(X,Y,Data);
    colormap(cmap)
    xlabel('x'); ylabel('y'); axis('equal'); 
    title(VarName);
    DivergingCentre( Data )
end


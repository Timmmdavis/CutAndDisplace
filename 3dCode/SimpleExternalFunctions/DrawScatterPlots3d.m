function DrawScatterPlots3d( X,Y,Z,cmap, varargin )
%DrawScatterPlots2d pass in X Y the colourmap and data and this plots each
%of the varargin for you as a 3D scatter plot. 
% X and Y are the X and Y Data.
% Cmap is the colourmap you will use
% I presume your variables are named well, the title takes the input arguments
% name

% Example: 
% DrawScatterPlots3d( X,Y,Z,cmap, S1,S2,S3 )

%   Copyright 2017, Tim Davis, The University of Aberdeen
for i=1:nargin-4
    VarName=inputname(i+4); %grabbing the input variables name, is a string
    Data=varargin{i};
    figure;
	scatter3(X(:),Y(:),Z(:),15,Data(:));
    colormap(cmap)
    xlabel('x'); ylabel('y'); axis('equal'); 
    title(VarName);
    DivergingCentre( Data )
end


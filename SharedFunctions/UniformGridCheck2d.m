function [ Uniform ] = UniformGridCheck2d( X,Y )
% UniformGridCheck2d: Checks X and Y grid points to see if these have a
%                   uniform linear spacing and they are shaped as such
%                   (i.e. each row/col increases by a uniform increment).
%                   Gridded data that has been transformed to col vects
%                   does therefore not fit this criteria. 
%
%                   Can be used to check before plotting with different
%                   functions or for using in interpolation. Returns 1 if
%                   true. Can use in if statements i.e.. if uniform spaced
%                   points plot with contourf else use scatter.
%   
% usage #1:
% [ Uniform ] = UniformGridCheck2d( X,Y )
%
% Arguments: (input)
% X,Y               - Data Points in X and Y. 
%
% Arguments: (output)
% Uniform           - Flag to say if the data points in X Y are spaced
%                    evenly on a grid ('1' means the data is uniform). 
%
% Example usage (1): Uniform grid
%
% [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
% [ Uniform ] = UniformGridCheck2d( X,Y )
% scatter(X(:),Y(:))
% 
% Example usage (2): Not uniform
%
% n=1000;
% xmv=0; ymv=0;
% x=rand(1,1000);
% y=rand(1,1000);
% X=x-n+xmv;
% Y=y-n+ymv;
% [X,Y]=ReshapeData2d( 100,10, X,Y )
% [ Uniform ] = UniformGridCheck2d( X,Y )
% scatter(X(:),Y(:))
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


if iscolumn(X)
    
    %If the data is col vecs we skip checking anything else
    flagx=0;
    flagy=0;
    
else    
    
    %Just checks gradient between points is 0.  Gradients where NaNs exist
    %are ignored.
    flagx = diff(X,[],1)==0 | isnan(diff(X,[],1));
    flagy = diff(Y,[],2)==0 | isnan(diff(Y,[],2));
 
end


Uniform= all(flagx(:)) && all(flagy(:)); %if the data is a uniform grid


end


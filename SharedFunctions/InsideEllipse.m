function [ External ] = InsideEllipse(x, y, rx, ry, tol_dist)
% InsideEllipse: Computes if points defined in vectors x/y are
%                   inside the ellipse specified by the semi-axis rx and
%                   ry. Returns a flag which is 1 when points are outside
%                   the ellipse.
%
% usage #1:
% [ External ] = InsideEllipse(x, y, rx, ry, tol_dist)
%
% Arguments: (input)
% X,Y               - The location of the point in 2D space. Can be matrix
%                    (grid) or vectors
%
% rx,ry             - The two semi-major axes aligned with the Cartesian
%                    plane.
%
% tol_dist          - The tolerance distance, adds some extra padding away
%                    from the ellipse.  
%
% Arguments: (output)
% External          - Flag to say if the point is inside (0) or outside (1)
%                    the ellipse.
%
% Example usage (1):
%
% x=1.1, y=1;
% rx=1; ry=1;
% tol_dist=0;
% [ External ] = InsideEllipse(x, y, rx, ry, tol_dist)
%
% Example usage (2):
% 
% [x,y]=meshgrid(-2:0.1:2,-2:0.1:2)
% rx=1; ry=1;
% tol_dist=0;
% [ External ] = InsideEllipse(x, y, rx, ry, tol_dist)
% scatter(x(:),y(:),[],External(:))
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Getting original size of X and Y (so we can return these to vectors
%afterwards).
Xorig=size(x);
%Turning into col vecs
x=x(:);y=y(:);
%Size of X. 
sz=numel(x);
%Our flag, assuming all are inside to start with
External = zeros(sz,1) ; 

%Doing for every point.
for i=1:sz
    if ((x(i)^2)/rx^2+(y(i)^2)/ry^2) > 1.0+tol_dist 
        External(i) = true ; 
    end  
end

%Reshaping back to the original size. 
External=reshape(External,Xorig);


end


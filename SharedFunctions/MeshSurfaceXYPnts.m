function [ Triangles,Points ] = MeshSurfaceXYPnts( X,Y )
%MeshSurfaceXYPnts: Given a set of points in plane X Y this meshes these to
%                   create a surface, we asssume these points lie close to
%                   the plane and are well distributed for meshing.
% usage #1:
% [ Triangles,Points ] = MeshSurfaceXYPnts( X,Y )
%
% Arguments: (input)
% X             - list of x points (vect)
%
% Y             - list of y points (vect)
%
% Arguments: (output)
% Points        - Columns 2 3 and 4 are the XYZ locations of one the
%                 corner points of a triangle. Column 1 is the index. 
%
% Triangles     - Triangles is a list where each row contains 3 index
%                 locations in "Points" which contains the XYZ location
%                 of each corner of the triangle.
% 
%
% Example usage:
%
% %Create a square of points
% x = [0 0 1 1]; 
% y = [0 1 1 0]; 
% %Mesh
% [ Triangles,Points ] = MeshSurfaceXYPnts( x,y );
% %Draw:
% PlotSlipDistribution3d(Triangles,Points,[],1 )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Getting rid of duplicate points
[ X,Y ] = RemoveDuplicatePoints2d( X,Y );

%Triangulating
TRI = delaunay(X,Y);
%Making points and triangles.
Triangles=TRI; clear TRI

%Appending points
Points=[(1:1:numel(X))',X(:),Y(:),zeros(size(X(:)))];


end


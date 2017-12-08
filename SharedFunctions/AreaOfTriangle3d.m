function [ Area,HPerim ] = AreaOfTriangle3d( x1,y1,z1,x2,y2,z2,x3,y3,z3 )
% AreaOfTriangle3d: Calculates the area of input triangles defined by the
%                   vectors of points that are corner 1,2 & 3 of the
%                   triangle. 
%
%                   Adapted from:
%                   %http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
%               
% usage #1:
% [ Area ] = AreaOfTriangle3d( x1,y1,z1,x2,y2,z2,x3,y3,z3 )
%
% Arguments: (input)
% x1,y1,z1          - X,Y & Z location of one corner of the triangle.
% x2,y2,z2          - X,Y & Z location of 2nd corner of the triangle.
% x3,y3,z3          - X,Y & Z location of 3rd corner of the triangle.
%
% Arguments: (output)
% Area              - The area of the triangle, will have the same input
%                    units as the X Y Z location of the points.
%
% HPerim            - Half the perimeter length of each triangle. 
%
% Example usage 1:
%
% %To find for a single flat lying right angle triangle:
% % y
% % ^  1
% % |  |\
% % |  | \
% % |  |  \
% % | 0 ¯ ¯ 1
% %  -------> x
%
% [ Area ] = AreaOfTriangle3d( 0,0,0,1,0,0,0,1,0 );
%
% Example usage 2:
%
% %To find for multiple triangles you would call as:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [P1,P2,P3] = CreateP1P2P3( Triangles,Points );
% [ Area,HPerim ] = AreaOfTriangle3d( P1(:,1),P1(:,2),P1(:,3),P2(:,1),P2(:,2),P2(:,3),P3(:,1),P3(:,2),P3(:,3) );
% PlotSlipDistribution3d(Triangles,Points,[],Area )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Calling internal function for segment distances (length between points or
%edge lengths of the triangles).
a=Distance3d(x1,y1,z1,x2,y2,z2);  
b=Distance3d(x2,y2,z2,x3,y3,z3); 
c=Distance3d(x3,y3,z3,x1,y1,z1) ;
%Now finding the area with internal function for Heron's formula.
[Area,HPerim] = HeronsFormula(a,b,c)  ;

function [Area,HPerim]=HeronsFormula(a,b,c)  
%Calculation for Heron's formula.    
%Inputs a b c are the edge lengths of the triangle    
    
%Half the length of the triangle perimeter.
HPerim=(a+b+c)./2;   
%Calculates the area.
Area=sqrt(HPerim.*(HPerim-a).*(HPerim-b).*(HPerim-c));    

end

function [Length]=Distance3d(x1,y1,z1,x2,y2,z2)   
%Internal function to find the edge lengths in 3D, Pythagoras therom.
%Inputs are xyz locations of two points in 3D space. 

%Squared Cartesian lengths.
SquaredLengths=(x1-x2).^2+(y1-y2).^2+(z1-z2).^2;
%Sqrt.
Length=sqrt(SquaredLengths); 

end

end


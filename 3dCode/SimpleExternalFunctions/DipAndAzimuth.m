function [Dip,Azimuth,Strike] = DipAndAzimuth( Triangles,Points,FaceNormalVector )
% DipAndAzimuth: Creates and plots the dip, azimuth and 
%                   strike of the surface in degrees. Plots the dip and
%                   azimuth of the surface and also exports these as column
%                   vectors. Result is independant of the surface's
%                   normals.
%               
% usage #1:
% [Dip,Azimuth,Strike] = DipAndAzimuth( Triangles,Points,FaceNormalVector )
%
% usage #2: (just drawing)
% DipAndAzimuth( Triangles,Points,FaceNormalVector )
%
% Arguments: (input)
% FaceNormalVector - The normal vector of the triangles of the surface.
%
%
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                     locations in "Points" which contains the XYZ location
%                     of each corner of the triangle.
%
% Arguments: (output)
%  Dip             - Dip of the surface (angle away from horizontal).
%
%  Azimuth         - Direction of the dip of the surface.
%
%  Strike          - Direction perpendicular to Azimuth (90 degrees counter
%                   clockwise).
%
% Example usage:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
% [Dip,Azimuth,Strike] = DipAndAzimuth( Triangles,Points,FaceNormalVector);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
CosAx=FaceNormalVector(:,1); 
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);

%Equations from Pollard Fletcher 2.102. Note these need some modification
%so they always give the dip and azimuth down dip.
Azimuth = atan2d(CosAx,CosAy); %extracting az (-180 to 180 degrees)
Dip = asind(-CosAz); %extracting dip, this is actually the dip of the normal

%Making az's 0-360;
Fullaz=Azimuth<0;
Azimuth(Fullaz)=Azimuth(Fullaz)+360;
%Now extracting positive and negative dips.
Convex=Dip<0; %Finding negative 'dip' values from pollard eq
Dip=90-abs(Dip); %Turning normal dip to triangle dip
%Flipping azimuth anywhere the dips where pointing up the way
Concave=Convex==0; %Creating a value to be used later
Azimuth(Concave,:)=Azimuth(Concave,:)+180;
NewVar=Azimuth>360;
Azimuth(NewVar)=Azimuth(NewVar)-360;

%Fault strike (degrees)
[ StrikeSlipCosine,~ ] = CalculateDSandSSDirs( FaceNormalVector );
Strike=atan2d(StrikeSlipCosine(:,1),StrikeSlipCosine(:,2));
FullStr=Strike<0;
Strike(FullStr)=Strike(FullStr)+360;


figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Azimuth);
xlabel('x'); ylabel('y'); axis('equal'); title('Azimuth ^{\circ}');colorbar;

figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Dip);
xlabel('x'); ylabel('y'); axis('equal'); title('Dip ^{\circ}');colorbar;

figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Strike);
xlabel('x'); ylabel('y'); axis('equal'); title('Strike ^{\circ}');colorbar;



end


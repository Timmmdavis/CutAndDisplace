function [Dip,Azimuth] = DipAndAzimuth( Triangles,Points,FaceNormalVector )
%Creates and plots the dip and azimuth of the surface
%   Plots the dip and azimuth of the surface and also exports these as
%   col vectors. The exported values always dip down, independant of
%   surface normals
%   Triangles is a list of the points that make the surface tris (col vec
%   N*3)
%   Points is a list of the point locations in XYZ and is a col vec (N*4)
%   where the first col is the index. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
CosAx=FaceNormalVector(:,1); 
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);

%Equations from Pollard Fletcher 2.102. Note these need some modification
%so they always give the dip and azimuth down dip.
AzimuthFromCosines = atan2d(CosAx,CosAy); %extracting az (-180 to 180 degrees)
DipFromCosines = asind(-CosAz); %extracting dip, this is actually the dip of the normal

%Making az's 0-360;
Fullaz=AzimuthFromCosines<0;
AzimuthFromCosines(Fullaz)=AzimuthFromCosines(Fullaz)+360;
%Now extracting positive and negative dips.
Convex=DipFromCosines<0; %Finding negative 'dip' values from pollard eq
DipFromCosines=90-abs(DipFromCosines); %Turning normal dip to triangle dip
%Flipping azimuth anywhere the dips where pointing up the way
Concave=Convex==0; %Creating a value to be used later
AzimuthFromCosines(Concave,:)=AzimuthFromCosines(Concave,:)+180;
NewVar=AzimuthFromCosines>360;
AzimuthFromCosines(NewVar)=AzimuthFromCosines(NewVar)-360;


Azimuth=AzimuthFromCosines;
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Azimuth);
xlabel('x'); ylabel('y'); axis('equal'); title('AzimuthFromCosines');colorbar;
Dip=DipFromCosines;
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Dip);
xlabel('x'); ylabel('y'); axis('equal'); title('DipFromCosines');colorbar;

%%%
%for fault strike you would call
% CosAx=FaceNormalVector(:,1);
% CosAy=FaceNormalVector(:,2);
% CosAz=FaceNormalVector(:,3);
% [ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector,CosAx,CosAy,CosAz );
% StrikeDegrees=acosd(StrikeSlipCosine(:,2));

end


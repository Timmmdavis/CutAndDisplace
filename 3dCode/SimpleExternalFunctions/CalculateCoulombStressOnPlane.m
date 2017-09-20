function [ CSS ] = CalculateCoulombStressOnPlane( X,Y,Z,FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Points,Triangles,cmap )
%CalculateCoulombStressOnPlane Calculates the coulomb stress change on an
%imported surface using the stress tensors at its midpoints. 
%Calculates the shear traction direction on the plane (every tri using a
%loop)
%Draws figures of tractions and stress.

%   Copyright 2017, Tim Davis, The University of Aberdeen

%X,Y,Z Midpoint locations
%FaceNormalVector = Direction cosines of plane
%Mu = Coeff Friction
%Sxx,Syy,Szz, Stress tensors
%Points = For drawing the surface
%Triangles= For drawing the surface
%cmap = For Coloring the CSS
%Cohesion = cohesive strength of fault surface

%Explict direction cosines
CosAx=FaceNormalVector(:,1);
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);

%Calculating the normal stresses on the second surface
[ Tnn ] = CalculateNormalTraction3d( Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz );

%Calculating traction on second surface 
[ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(  Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz );

% %Total traction vector
T=sqrt((Tx.^2)+(Ty.^2)+(Tz.^2)); %TractionVectorOnPlane

%Pollard and fletcher Book Eq 6.53 %Max Shear stress
Ts_maxShr=sqrt((abs(T).^2)-(abs(Tnn).^2));

%CoulombShearStress %Pollard and fletcher Book Eq  9.40 
CSS=abs(Ts_maxShr)+(Mu.*Tnn)-Cohesion;

%Creating and drawing Shearing traction vector
TsVector=zeros(size(FaceNormalVector));

%Equations below taken in part from
%Professor Ramón Arrowsmith's online lectures
%"arrowsmith510.asu.edu/TheLectures/Lecture16/Lecture16_3Dstress.ppt"
%Now for the shear traction; use the McKenzie construction
for i=1:size(FaceNormalVector(:,1))
N=FaceNormalVector(i,:)';
T=[Tx(i,:);Ty(i,:);Tz(i,:)];
B = cross(T,N);         %vector normal to the plane containing T and N
Ts = cross(N,B);        %shear traction direction
Ts_mag=Ts_maxShr(i);        %Shear traction magnitude
Ts(1) = Ts(1)./Ts_mag;  %X dir
Ts(2) = Ts(2)./Ts_mag;  %Y dir
Ts(3) = Ts(3)./Ts_mag;  %Z dir
TsVector(i,:)=Ts';
end

%Drawing total traction and normal stress
figure;quiver3(X,Y,Z,Tx(:,1),Ty(:,1),Tz(:,1))
xlabel('x'); ylabel('y'); axis('equal'); title('Total Traction Vector');
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Tnn);colormap(cmap)
DivergingCentre( Tnn )
hold off

%Drawing shear traction and CSS 
figure;quiver3(X,Y,Z,TsVector(:,1),TsVector(:,2),TsVector(:,3))
xlabel('x'); ylabel('y'); axis('equal'); title('Shear Traction Vector and CSS');
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),CSS);colormap(cmap)
DivergingCentre( CSS )
hold off



end


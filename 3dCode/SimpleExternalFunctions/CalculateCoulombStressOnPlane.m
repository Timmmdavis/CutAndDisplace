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

%Calculates the directions of the dipslip and ss directions
[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector,CosAx,CosAy,CosAz );

%TractionStrikeSlip
[ Tss ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,StrikeSlipCosine );
%TractionStrikeSlip
[ Tds ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine );

%Cart vector components of the strike slip traction 
TssCart=bsxfun(@times,StrikeSlipCosine,Tss);
%Cart vector components of the dip slip traction 
TdsCart=bsxfun(@times,DipSlipCosine,Tds);

%3D vector addition of these Cart components (normalised as we want the
%vector). This is the max shear vector direction.
TsVector=normr(TssCart+TdsCart);

%Total traction vector on the plane
T=sqrt((Tx.^2)+(Ty.^2)+(Tz.^2)); 

%Max Shear stress %Pollard and Fletcher Book Eq 6.53 
Ts_maxShr=sqrt((abs(T).^2)-(abs(Tnn).^2));

%Function to calculate CSC
[ CSS ] = CalculateCoulombShearStress( Tnn,Ts_maxShr,Mu,Cohesion );

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


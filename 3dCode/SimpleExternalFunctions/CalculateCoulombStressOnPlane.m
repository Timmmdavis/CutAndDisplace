function [ CSS,TsMaxShr,TsMaxShrDir ] = CalculateCoulombStressOnPlane(MidPoint,FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Points,Triangles,cmap )
% CalculateCoulombStressOnPlane: Calculates the Coulomb stress change on an
%                   imported triangulated surface using the stress tensors
%                   at its midpoints. Calculates the shear traction
%                   direction on the plane and draws figures of tractions
%                   and stress.
%               
% usage #1:
%[ CSS ] = CalculateCoulombStressOnPlane...
%( MidPointObs,FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Points,Triangles,cmap )
%
% Arguments: (input)
%   MidPoint       - A 3*n vector that is the XYZ locations of the
%                   midpoints of the imported fault surface that you find
%                   the Coulomb stress change on.
%
% FaceNormalVector - The normal vector of the triangles of the surface.
%
% Sxx,Syy,Szz...
% Sxy,Sxz,Syz      - The calculated stress tensor at each midpoint of the
%                   surface.
%
% Mu                - The coefficient of friction a every midpoint of the
%                    surface. (column vec)
%
% Cohesion          - The Cohesive strength of the surface (or sliding
%                    friction).
%
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                     locations in "Points" which contains the XYZ location
%                     of each corner of the triangle.
%
% cmap              -  A colourmap that MATLAB can use. See func
%                   "colormap_cpt.m" to produce one. 
%
% Arguments: (output)
%       CSS        - The Coulomb stress change at each point
%
%  TsMaxShr        - The magnitude of the maxiumum shear stress at each
%                   point. 
%
%  TsMaxShrDir     - Direction cosines of the maximum shear direction
%                   (CosAx,CosAy,CosAz)
%
% Example usage:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
% Sxx=ones(size(FaceNormalVector(:,1)));
% Syy=zeros(size(FaceNormalVector(:,1)));
% Szz=Syy;
% Sxy=Syy;
% Sxz=Syy;
% Syz=Syy;
% Mu=0.6; 
% Cohesion=0.2;
% [ CSS,TsMaxShr,TsMaxShrDir ] = CalculateCoulombStressOnPlane...
% (MidPoint,FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Points,Triangles,[]);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Grabbing the midpoints of the surface
X=MidPoint(:,1);
Y=MidPoint(:,2);
Z=MidPoint(:,3);

%Explict direction cosines
CosAx=FaceNormalVector(:,1);
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);

%Calculating the normal stresses on the second surface
[ Tn ] = CalculateNormalTraction3d( Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz );

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
TsMaxShrDir=normr(TssCart+TdsCart);

%Total traction vector on the plane
T=sqrt((Tx.^2)+(Ty.^2)+(Tz.^2)); 

%Max Shear stress %Pollard and Fletcher Book Eq 6.53 
TsMaxShr=sqrt((abs(T).^2)-(abs(Tn).^2));

%Function to calculate CSC
[ CSS ] = CalculateCoulombShearStress( Tn,TsMaxShr,Mu,Cohesion );

%Drawing total traction and normal stress
figure;quiver3(X,Y,Z,Tx(:,1),Ty(:,1),Tz(:,1))
xlabel('x'); ylabel('y'); axis('equal'); title('Total Traction Vector');
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Tn);
if isempty(cmap); colormap('default'); else; colormap(cmap); end %draw with the imported value
DivergingCentre( Tn )
hold off

%Drawing shear traction and CSS 
figure;quiver3(X,Y,Z,TsMaxShrDir(:,1),TsMaxShrDir(:,2),TsMaxShrDir(:,3))
xlabel('x'); ylabel('y'); axis('equal'); title('Shear Traction Vector and CSS');
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),CSS);
if isempty(cmap); colormap('default'); else; colormap(cmap); end %draw with the imported value
DivergingCentre( CSS )
hold off

end


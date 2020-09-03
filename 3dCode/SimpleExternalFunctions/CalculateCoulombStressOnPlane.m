function [ CSS ] = CalculateCoulombStressOnPlane(FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Rake )
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
% Rake              - Defined as clockwise away from the pure thrust direction
%                     facing down the normal onto the fault. Supplied in
%                     degrees.
%
% Arguments: (output)
%       CSS        - The Coulomb stress change at each point
%
% Example usage A) on a mesh:
%
% Rake=0; %Thrust
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles,0);
% Sxx=zeros(size(FaceNormalVector(:,1)));
% Syy=zeros(size(FaceNormalVector(:,1)));
% Szz=ones(size(FaceNormalVector(:,1)));
% Sxy=Sxx;
% Sxz=Sxx;
% Syz=Sxx;
% Mu=0.6; 
% Cohesion=0.2;
% [ CSS ] = CalculateCoulombStressOnPlane...
% (FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Rake);
% hold on
% trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),CSS,'LineStyle','none'); WhiteFigure;
% DivergingCentre( CSS ); axis('equal'); colormap('cool');
% hold off
%
% Example usage B) on a series of points with imaginary planes:
%
% %Creating grid with user defined sampling
% X = linspace(-10,15,25); 
% Y = linspace(-10,15,25); 
% [X,Y] = meshgrid(X,Y); 
% [dimx,dimy] = size(X);  
% %Planes at each point
% SurfaceDip=89;
% SurfaceAzimuth=0;
% Rake=270; %Thrust
% %Normal vector facing upwards (flat plane)
% FaceNormalVectorObs=[0, 0, 1];
% %First we rotate around dip
% [FaceNormalVectorObs]=RotateCosine3d(FaceNormalVectorObs,deg2rad(-SurfaceDip),'x');
% %now around azimuth
% [FaceNormalVectorObs]=RotateCosine3d(FaceNormalVectorObs,deg2rad(SurfaceAzimuth),'z');
% %Repeat for every point (could have different if needed) 
% FaceNormalVectorObs=repmat(FaceNormalVectorObs,numel(X),1);
% Sxx=zeros(size(FaceNormalVectorObs(:,1)));
% Syy=zeros(size(FaceNormalVectorObs(:,1)));
% Szz=zeros(size(FaceNormalVectorObs(:,1)));
% Sxy=ones(size(FaceNormalVectorObs(:,1)));
% Sxz=Syy;
% Syz=Syy;
% Mu=0.6; 
% Cohesion=0.2;
% [ CSS ] = CalculateCoulombStressOnPlane...
% (FaceNormalVectorObs,Sxx,Syy,Szz,Sxy,Sxz,Syz,Mu,Cohesion,Rake);
% DrawContourFPlots2d( X,Y,[],reshape(CSS,size(X)) );
%   
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

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
%TractionDipSlip
[ Tds ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine );

Ts=Tds*cos(deg2rad(Rake))+Tss*sin(deg2rad(Rake));

%Function to calculate CSC
[ CSS ] = CalculateCoulombShearStress( Tn,Ts,Mu,Cohesion );
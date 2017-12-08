function [Ux,Uy,Uz]=Mogi1958_SphericalCavity(D,X,Y,Radius,P,nu,mu)
% Mogi1958_SphericalCavity: Returns the displacement of the ground surface
%               (Cartesian) due to a pressurised spherical source at depth.
%               The ground surface is approximately traction free.
%
%               See: Segall, P., 2010. Earthquake and volcano deformation.
%               Princeton University Press. & Lisowski, M., 2007.
%               Analytical volcano deformation source models. In Volcano
%               Deformation (pp. 279-304). Springer Berlin Heidelberg.
%
% usage #1:
% [Ux,Uy,Uz]=Mogi1958_SphericalCavity(D,X,Y,Radius,P,nu,mu)
%
% Arguments: (input)
%       Radius - The cavities radius
%
%       D      - The depth of the centre of the cavity.
%
%       X,Y    - The location of the observation points at the ground.
%               surface (in relation to the chambers centre). 
%
%       P      - The internal pressure within the hole. 
%
%       nu    - Poisson's Ratio.
%
%       mu    - Shear Modulus.
%
%
% Arguments: (output)
% Ux,Uy,Uz    - The Cartesian components of displacement at the ground
%              surface.
%
% Example usage:
%
%  [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  a = 1;
%  P = 1; 
%  D=10; 
%  Radius=1; 
%  nu=0.25;
%  mu=50;
%  
%  [Ux,Uy,Uz]=Mogi1958_SphericalCavity(D,X,Y,Radius,P,nu,mu);
% 
%  scl=150000;
%  mesh(X+(Ux.*scl),Y+Uy.*scl,zeros(size(X))+Uz.*scl);
%  hold on
%  [x,y,z] = ellipsoid(0,0,0,Radius,Radius,Radius,20);
%  hSurface=surf(x, y, z-D);
%  axis('equal');
%  caxis([min(Uz(:).*scl),max(Uz(:).*scl)])
%  title('Exaggerated ground deformation over mogi point source')
% 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Distance from centre of chamber to point.
R=sqrt(D.^2+X.^2+Y.^2); 

%Right side of eq. 8.15 (Lisowski Book)
lsi=Radius^3*P*((1-nu)/mu);

%Cartesian components of displacement
Ux=(X./(R.^3)).*lsi;
Uy=(Y./(R.^3)).*lsi;
Uz=(D./(R.^3)).*lsi;



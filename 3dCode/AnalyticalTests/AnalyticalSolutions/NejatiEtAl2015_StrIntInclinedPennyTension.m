function [K1,K2,K3] = NejatiEtAl2015_StrIntInclinedPennyTension(Szz,Beta,Theta,Radius,nu)
% NejatiEtAl2015_StrIntInclinedPennyTension: Returns the stress intensity 
%               factors at the tips of a penny shaped crack under tension.
%				Equation D1 in: Nejati, M., Paluszny, A. and Zimmerman,
%				R.W., 2015. A disk-shaped domain integral method for the
%				computation of stress intensity factors using tetrahedral
%				meshes. International Journal of Solids and Structures, 69,
%				pp.230-251.
%
% usage #1: 
% [K1,K2,K3] = NejatiEtAl2015_StrIntInclinedPennyTension(Szz,Beta,Radius)
%
% Arguments: (input)
% Szz             - The remote uniaxial tension along the Z-axis.
%
% Beta            - The inclination of the crack (relative to Z-axis). Note
%                   this is inclined in the XZ-axis.
%
% Theta           - The location on the crack wall, measured anticlockwise 
%                   away from X-axis (in crack plane), degrees. 
%
% Radius          - The radius of the crack.
%
% nu              - Poisson's ratio of material
%
% Arguments: (output)
% K1              - K1 at the crack tips.
%
% K2              - K2 at the crack tips.
%
% K3              - K3 at the crack tips.
%
% Example usage:
%
% %uniaxial tension in z
% Szz=1; 
% Radius=1; 
% [ x,y ] = CreateCircleXY( 100,Radius );
% z=zeros(size(x));
% Theta=cart2pol(x,y)
% %Rotate this (XZ)
% Beta=20;
% [x,z] = RotateObject2d(x,z,deg2rad(90-Beta));
% nu=0.25;
% [K1,K2,K3] = NejatiEtAl2015_StrIntInclinedPennyTension(Szz,Beta,Theta,Radius,nu)
% scatter3(x,y,z,[],K2);
% axis('equal')
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University

% %%Inclinded Penny crack Eq D1. 
% %Nejati M, Paluszny A, Zimmerman RW, 2015
K1=ones(size(Theta))*(2*Szz*sqrt(Radius/pi)*sind(Beta)^2);
% K2=((2*Szz*sqrt(Radius/pi))/(2-nu))*sind(2*Beta)*cosd(Theta);
% K3=((2*(1-nu)*Szz*sqrt(Radius/pi))/(2-nu))*sind(2*Beta)*sind(Theta);

%Penny crack tension Tada Stress analysis handbook P.24.22:
K2=(4/(pi*(2-nu)))*(Szz*sind(Beta)*cosd(Beta))*cosd(Theta)*sqrt(pi*Radius);
K3=((4*(1-nu))/(pi*(2-nu)))*(Szz*sind(Beta)*cosd(Beta))*sind(Theta)*sqrt(pi*Radius);


% % If the penny crack is flat this reduces to:
% K1=2*Szz*sqrt(Radius/pi);




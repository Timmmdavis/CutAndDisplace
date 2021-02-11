function [Ux,Uy] = K1DispDistro(K1, X,Y, TipDip,mu,nu)
%K1StressDistro Distribution of displacements around the tip of a crack due to
%stress intensity K1
%
% usage #1:
% [Ux,Uy] = K1DispDistro(K1, X,Y, TipDip,mu,nu)
%
% Arguments: (input)
%  K1   - Value of K1 at the cracks tip (just a float not a matrix)
%
%  X    - X distance of point from the tip of the crack (can be a matrix, vector or float)
%
%  Y    - Y distance from the tip of the crack (can be a matrix, vector or float)
%
% TipDip - Dip angle of the tip in degrees
%
% mu    - Shear modulus
%
% nu    - Poisson's ratio
%
% Arguments: (output)
% Ux,Uy - Displacements tensors at the points around the crack, back in the
% global coordinate system
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%[dimx,dimy] = size(X);  

%Rotate so we are in the tip coordiantes:
[Xrot,Yrot] = RotateObject2d(X(:),Y(:),deg2rad(-TipDip));
[theta,r]=cart2pol(Xrot,Yrot);

%Page 3 Tada Stress analysis of cracks handbook
Cons1=K1./mu;
Cons2=sqrt(r./(2*pi));
Ux=Cons1.*Cons2.*(cos(theta./2).*(1-(2*nu)+sin(theta./2).^2));
Uy=Cons1.*Cons2.*(sin(theta./2).*(2-(2*nu)-cos(theta./2).^2));

%Convert back to global coords
[Ux,Uy] = RotateObject2d(Ux(:),Uy(:),deg2rad(TipDip));

end


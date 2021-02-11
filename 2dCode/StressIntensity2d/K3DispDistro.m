function [Uz] = K3DispDistro(K3, X,Y, TipDip,mu,nu)
%K3StressDistro Distribution of displacements around the tip of a crack due to
%stress intensity K3
%
% usage #1:
% [Uz] = K3DispDistro(K3, X,Y, TipDip,mu,nu)
%
% Arguments: (input)
%  K3   - Value of K3 at the cracks tip (just a float not a matrix)
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
Cons1=K3./mu;
Cons2=sqrt(r./(2*pi));
Uz=Cons1.*Cons2.*(sin(theta./2));

end


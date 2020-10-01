function [Sxz,Syz] = K3StressDistro(K3, X,Y, TipDip)
%K1StressDistro Distribution of stresses around the tip of a crack due to
%stress intensity K1
%
% usage #1:
% [Sxx,Syy,Sxy] = K2StressDistro(K2, X,Y, TipDip)
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
% Arguments: (output)
% Sxz,Syz - Stress tensors at the points around the crack, back in the
% global coordinate system. These are the out of plane components.
%
% Example usage:
% %Reproduce P&F Fig9.31 box c)
% K3=1
% theta=linspace(-pi,pi,100);
% r=1;
% [X,Y] = pol2cart(theta,r);
% [Sxz,Syz] = K3StressDistro(K3, X,Y, 0);
% figure;hold on
% plot(theta,Sxz)
% plot(theta,Syz)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%[dimx,dimy] = size(X);  

%Rotate so we are in the tip coordiantes:
[Xrot,Yrot] = RotateObject2d(X(:),Y(:),deg2rad(-TipDip));
[theta,r]=cart2pol(Xrot,Yrot);

%Eq 9.72 Pollard and Fletcher
Cons1=K3./sqrt(2*pi*r);
Sxz=Cons1.*(-sin(theta/2));
Syz=Cons1.*(cos(theta/2));

%Convert back to global coords - No need!
%CosAx=cos(TipDip);
%CosAy=sin(TipDip);
%[ Sxz,Syz ] = StressTensorTransformation2d(Sxz,Syz,zeros(size,CosAx,CosAy );

%Back to input grid size
%[Sxx,Syy,Sxy]=ReshapeData2d...
%    ( dimx,dimy,Sxx,Syy,Sxy );

end


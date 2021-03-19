function [Sxx,Syy,Sxy] = K1StressDistro(K1, X,Y, TipDip)
%K1StressDistro Distribution of stresses around the tip of a crack due to
%stress intensity K1
%
% usage #1:
% [Sxx,Syy,Sxy] = K1StressDistro(K1, X,Y, TipDip)
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
% Sxx,Syy,Sxy - Stress tensors at the points around the crack, back in the
% global coordinate system
%
% Example usage:
% %Reproduce P&F Fig9.31 box a)
% K1=1
% theta=linspace(-pi,pi,100);
% r=1;
% [X,Y] = pol2cart(theta,r);
% [Sxx,Syy,Sxy] = K1StressDistro(K1, X,Y, 0);
% figure;hold on
% plot(theta,Sxx)
% plot(theta,Syy)
% plot(theta,Sxy)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%[dimx,dimy] = size(X);  

%Rotate so we are in the tip coordiantes:
[Xrot,Yrot] = RotateObject2d(X(:),Y(:),deg2rad(-TipDip));
[theta,r]=cart2pol(Xrot,Yrot);

%Eq 9.72 Pollard and Fletcher
Cons1=K1./sqrt(2*pi*r);
Sxx=Cons1.*(cos(theta/2).*(1-sin(theta/2).*sin((3*theta)/2)));
Syy=Cons1.*(cos(theta/2).*(1+sin(theta/2).*sin((3*theta)/2)));
Sxy=Cons1.*(sin(theta/2).*cos(theta/2).*cos((3*theta)/2));

%Convert back to global coords
CosAx=cos(deg2rad(TipDip));
CosAy=sin(deg2rad(TipDip));
[ Sxx,Syy,Sxy ] = StressTensorTransformation2d(Sxx,Syy,Sxy,CosAx,CosAy );

%Back to input grid size
%[Sxx,Syy,Sxy]=ReshapeData2d...
%    ( dimx,dimy,Sxx,Syy,Sxy );

end


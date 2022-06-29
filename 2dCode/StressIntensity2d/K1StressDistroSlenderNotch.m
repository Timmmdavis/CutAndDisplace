function [Sxx,Syy,Sxy] = K1StressDistroSlenderNotch(K1, X,Y, TipDip, rho)
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
% rho - Radius of curvature at tip...
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

Cons1=K1./sqrt(2*pi*r);
Cons2=rho./(2.*r);

%Eq.4 Tada 2000 - first term
Sxxterm1=Cons1.*Cons2.*(-cos((3*theta)/2));
Syyterm1=Cons1.*Cons2.*(cos((3*theta)/2));
Sxyterm1=Cons1.*Cons2.*(-sin((3*theta)/2));

%Eq 9.72 Pollard and Fletcher
Sxxterm2=Cons1.*(cos(theta/2).*(1-sin(theta/2).*sin((3*theta)/2)));
Syyterm2=Cons1.*(cos(theta/2).*(1+sin(theta/2).*sin((3*theta)/2)));
Sxyterm2=Cons1.*(sin(theta/2).*cos(theta/2).*cos((3*theta)/2));

%Sum terms
Sxx=Sxxterm1+Sxxterm2;
Syy=Syyterm1+Syyterm2;
Sxy=Sxyterm1+Sxyterm2;


%Convert back to global coords
CosAx=cos(deg2rad(TipDip));
CosAy=sin(deg2rad(TipDip));
[ Sxx,Syy,Sxy ] = StressTensorTransformation2d(Sxx,Syy,Sxy,CosAx,CosAy );

%Back to input grid size
%[Sxx,Syy,Sxy]=ReshapeData2d...
%    ( dimx,dimy,Sxx,Syy,Sxy );

end


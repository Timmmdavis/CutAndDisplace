function [Srr,Stt,Srt] = K2StressDistroRadialCoords(K2, X,Y, TipDip)
%K2StressDistro Distribution of stresses around the tip of a crack due to
%stress intensity K2
%
% usage #1:
% [Srr,Stt,Srt] = K2StressDistro(K2, X,Y, TipDip)
%
% Arguments: (input)
%  K2   - Value of K2 at the cracks tip (just a float not a matrix)
%
%  X    - X distance of point from the tip of the crack (can be a matrix, vector or float)
%
%  Y    - Y distance from the tip of the crack (can be a matrix, vector or float)
%
% TipDip - Dip angle of the tip in degrees
%
% Arguments: (output)
% Srr,Stt,Srt - Radial stress tensors at the points around the crack
%
% Example usage:
% K2=1
% theta=linspace(-pi,pi,100);
% r=1;
% [X,Y] = pol2cart(theta,r);
% [Srr,Stt,Srt] = K2StressDistro(K2, X,Y, 0);
% figure;hold on
% plot(theta,Srr
% plot(theta,Stt)
% plot(theta,Srt)
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%[dimx,dimy] = size(X);  

%Rotate so we are in the tip coordinates:
[Xrot,Yrot] = RotateObject2d(X(:),Y(:),deg2rad(-TipDip));
[theta,r]=cart2pol(Xrot,Yrot);

%Eq 9.72 Pollard and Fletcher
Cons1=K2./sqrt(2*pi*r);
Srr=Cons1.*(sin(theta/2).*(1-3*sin(theta/2).^2));
Stt=Cons1.*(-sin(theta/2).*(3*cos(theta/2).^2));
Srt=Cons1.*(cos(theta/2).*(1-3*sin(theta/2).^2));


%Back to input grid size
%[Sxx,Syy,Sxy]=ReshapeData2d...
%    ( dimx,dimy,Sxx,Syy,Sxy );

end


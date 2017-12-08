function [Sxx,Syy,Sxy,Srr,Stt,Srt]=Kirsch1898_PressurisedHole(X,Y,a,P,Sxx,Syy)
% Kirsch1898_PressurisedHole: Returns the stress tensors components around
%               a hole pressurised for peturbed by a remote stress.
%
%               See: Pollard, D.D. and Fletcher, R.C., 2005. Fundamentals
%               of structural geology. Cambridge University Press.
%
% usage #1:
% [Sxx,Syy,Sxy,Srr,Stt,Srt]=Kirsch1898_PressurisedHole(X,Y,a,P,Sxx,Syy)
%
% Arguments: (input)
%       a      - The central holes radius
%
%       X,Y    - The location of the observation points surronding the
%               hole.
%
%       P      - The internal pressure within the hole. 
%
%     Sxx,Syy  - The remote stresses at infinity.  
%
%
% Arguments: (output)
% Sxx,Syy,Sxy  - The Cartesian components of stress at points X,Y
%               surronding the hole
%
% Srr,Stt,Srt  - The radial components stress at points X,Y
%               surronding the hole
%
%
% Example usage:
%
%
%  [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  a = 1;
%  P = 1; 
%  Sxx = -1; 
%  Syy = -1;
%  
%  [Sxx,Syy,Sxy,Srr,Stt,Srt]=Kirsch1898_PressurisedHole(X,Y,a,P,Sxx,Syy);
% 
%  DrawContourFPlots2d( X,Y,[], Sxx,Syy,Sxy,Srr,Stt,Srt)
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
%
% Script modified from from Pollard/Fletcher Fundamentals of Structural Geology
% Don't reproduce without the above statement ^^^

if nargin<5
    Sxx=0;
    Syy=0;
end    

%Internal pressure in hole
P=P(1,1);

%Radial coordiantes for the input observation point locations
[theta,rho] = cart2pol(X,Y);

%Creating some common variables
S2t = sin(2*theta); 
C2t = cos(2*theta); 
r2 = (a./rho).^2; 
r4 = r2.^2;

% Polar stress components
Srr = -(0.5*(-Sxx-Syy)*(1-r2))+(P*r2)-(0.5*(-Sxx+Syy)*((1-4*r2+3*r4).*C2t));
Stt = -(0.5*(-Sxx-Syy)*(1+r2))-(P*r2)+(0.5*(-Sxx+Syy)*((1+3*r4).*C2t));
Srt = 0.5*(-Sxx+Syy)*((1+2*r2-3*r4).*S2t);

%Removing points inside of the hole 
Tol=0;
Srr(rho<a-Tol) = nan; Stt(rho<a-Tol) = nan; Srt(rho<a-Tol) = nan;

%Converting the components into Cartesian tensor
dimx = size(Srr,1);
dimy = size(Srr,2);
%Calling external function
[Sxx,Syy,Sxy] = StressTensorTransformation2d(Srr(:),Stt(:),Srt(:),cos(theta(:)),sin(theta(:)));
[Sxx,Syy,Sxy] = ReshapeData2d(dimx,dimy,Sxx,Syy,Sxy);


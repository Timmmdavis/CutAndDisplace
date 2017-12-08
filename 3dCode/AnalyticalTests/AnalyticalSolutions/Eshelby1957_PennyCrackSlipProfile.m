function [Disp]=Eshelby1957_PennyCrackSlipProfile(mu,nu,Tn,Ts,r)
% Eshelby1957_PennyCrackSlipProfile: Returns the displacement
%               of fracture walls for a penny shaped crack
%               loaded by a constant shear/normal traction. Crack half
%               length here is 1.
%               
%               See: Segall, P., 2010. Earthquake and volcano deformation. 
%               Princeton University Press.
%
% usage #1:
% [Disp]=Eshelby1957_PennyCrackSlipProfile(mu,nu,Tn,Ts,r)
% usage #2: (Assuming this is lies in the XY plane
% [Disp]=Eshelby1957_PennyCrackSlipProfile(mu,nu,Syy,Sxz/Syz,r)
%
% Arguments: (input)
%       nu    - Poisson's Ratio.
%
%       mu    - Shear Modulus.
%
%       Tn    - Normal traction on the fractures face. 
%
%       Ts    - Shear traction on the fractures face. 
%
%       r     - The location of the points on the fractures walls.(radial
%              distance from centre).Assuming a=1 these should be from -1
%              to 1.
%
% Arguments: (output)
%       Disp   - The Burgers vector (separation of dislocation walls) at
%               each discrete r location.
%
% Example usage:
%
%  r  = linspace(-1,1,50); 
%  mu = 500; 
%  nu = 0.25;
%  Ts= 1; 
%  Tn= 1;
% 
%  [Disp]=Eshelby1957_PennyCrackSlipProfile(mu,nu,Tn,Ts,r)
% 
%  plot(r,Disp); xlabel('radial-coord'); ylabel('Displacement of walls') 
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


% Unit half-length displacement discontinuity. This is going to have its 
% centre at 0,0. 
a = 1;  

% Calculate displacement of penny shaped crack walls.
% These are seperated into seperate components based on related driving
% stress so the equation as a whole is easier to read.

%Shearing of crack walls % %One side of crack (eq 4.74 segall)
us=((4*(1-nu)*a*Ts)/(pi*(2-nu)*mu))*sqrt(1-((r.^2)/(a^2)));
%Opening of crack walls % %Normal displacement, (eq 4.73 segall)
ut=((2*a*(1-nu)*Tn)/(pi*mu))*sqrt(1-((r.^2)/(a^2)));
%Both (pythag therom)
Disp=sqrt((us.^2)+(ut.^2));

% distance between relative faces;
Disp=Disp*2;


% %%%Pressure to vol relationship for crack loaded by shear stress
% Vol=(8*Ts*a^3*(nu - 1))/(3*mu*(nu - 2))
% %%%Pressure to vol relationship for crack loaded by normal stress
% Vol=(2*Tn*a^3*(1 - r^2/a^2)^(1/2)*(nu - 1))/mu
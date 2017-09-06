function [Disp]=Func_SlipPollardSegall(G,PR,Syy,Sxy,x)

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Assuming the crack lies in the X axis and only 1 stress acts on the crack
%at a time
%G shear mod
%PR nu (Poisson's ratio)
%Syy driving stress normal to the crack
%Sxy driving stress shearing the crack
%x location of points on the crack wall. Assuming a=1 these should be from
%-1 to 1

% Unit half-length displacement discontinuity. This is going to have its 
% centre at 0,0. 
a = 1;  

% Calculate displacement eq 8.35a and 8.35b Pollard and Segall 1984
% These are seperated into seperate components based on related driving stress so the equation as a whole is easier to read 
Disp=(2*(1-PR)/G)*sqrt(a.^2-x.^2)*(Syy+Sxy);
Mode3Disp=((2*(1-PR)/G)*sqrt(a.^2-x.^2))/1-PR;


%%%Pressure vol relationship
% %Elastic constant
% Vol=-(pi*(Tn+Ts)*a^2*(nu - 1))/mu
function [Disp]=Func_SlipPennyCrackSneddon(mu,nu,Tn,Ts,r)

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Assuming the crack lies in the X axis and only 1 stress acts on the crack
%at a time
%mu shear mod
%nu (Poisson's ratio)
%Tn driving stress normal to the crack
%Ts driving stress shearing the crack
%r location of points on the crack wall (radial). Assuming a=1 these should be from
%-1 to 1

% Unit half-length displacement discontinuity. This is going to have its 
% centre at 0,0. 
a = 1;  

% Calculate displacement of penny shaped crack walls
% These are seperated into seperate components based on related driving stress so the equation as a whole is easier to read 

%Shearing of crack walls % %One side of crack (eq 4.74 segall)
us=((4*(1-nu)*a*Ts)/(pi*(2-nu)*mu))*sqrt(1-((r.^2)/(a^2)));
%Opening of crack walls % %Normal displacement, (eq 4.73 segall)
ut=((2*a*(1-nu)*Tn)/(pi*mu))*sqrt(1-((r.^2)/(a^2)));
%Both (pythag therom)
Disp=sqrt((us.^2)+(ut.^2));

% distance between relative faces;
Disp=Disp*2;


% %%%Pressure to vol relationship for crack loaded by normal stress
% VPredicted=(((Tn*4)/mu)+((-5-1/3)*nu+(1+1/3)))*(a^3);
% TnPred=(((VPredicted./(a^3))-((-5-1/3)*nu+(1+1/3)))*mu)/4;
% test=1;
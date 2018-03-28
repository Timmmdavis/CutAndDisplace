function [K1,K2] = ShiahLin2004_StrIntCurvedCrack(Pressure,CrkHeight,HalfLength)
% ShiahLin2004_StrIntCurvedCrack: Returns the stress intensity factor at
%               the tips of a curved crack with internal pressure (or
%               subject to internal pressure).
%				Equation 34 of: Shiah, Y.C. and Lin, Y.J., 2004. An
%				infinitely large plate weakened by a circular-arc crack
%				subjected to partially distributed loads. Journal of 
%               engineering mathematics, 49(1), pp.1-18.
% 
%
% usage #1: 
% [K1,K2] = ShiahLin2004_StrIntCurvedCrack(Pressure,CrkHeight,HalfLength)
%
% Arguments: (input)
% Pressure        - The pressure inside the crack (also corresponds to the
%                   remote stresses (biaxial tension). 
%
% CrkHeight       - Radius minus the length from the centre of the problem 
%                   to the cracks maximum height ('?L' in figure 2 of the 
%                   paper). Really the length the crack spans in Y. 
%
% HalfLength      - Half the width of the crack in the X-axis. 
%
% Arguments: (output)
% K1              - K1 at the crack tips
%
% K2              - K2 at the crack tips
%
% Example usage:
%
% P=1; %Pressure (or biaxial tension)
% r=1; %(L in sih)
% HalfLength=0.5;
% x=linspace(-HalfLength,HalfLength,50);    
% y=-(sqrt(r.^2-(x).^2));
% plot(x,y,'r');hold on;scatter(0,0,'k');
% axis('equal')
% %Calculate span of crack in Y. 
% CrkHeight=max(y)-min(y);
% [K1,K2] = ShiahLin2004_StrIntCurvedCrack(P,CrkHeight,HalfLength)
% text(max(x),max(y),strcat({'K1=',num2str(K1),'K2=',num2str(K2)}))
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University

K1=sqrt(pi*HalfLength)*sqrt(1-(CrkHeight.^2))/(1+(CrkHeight.^2))*Pressure;
K2=sqrt(pi*HalfLength)*(CrkHeight/(1+(CrkHeight.^2)))*Pressure;

% % If crack is flat this reduces to:
% K1=Pressure*sqrt(pi*HalfLength);




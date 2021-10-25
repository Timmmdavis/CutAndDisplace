function [K1,K2] = StressIntensity2d(EndElsLoc,Dn,Ds,HalfLength,E,nu...
    ,LineNormalVector,MidPoint,P1,P2)
%StressIntensity2D: Returns the value of stress intensity at the end elements.
%               Outputs are vectors the size of Dn and Dn with K1 and K2 at
%               locations described in the EndElsLoc flag and zeros
%               everywhere else.
%               
% usage:
% [K1,K2] = StressIntensity2d(EndElsLoc,Dn,Ds,HalfLength,E,nu...
%    ,LineNormalVector,MidPoint,P1,P2)
%
% Arguments: (input)
% EndElsLoc          - A vector the same length Dn etc. Acts as a flag to
%                      say if an element is an endpoint. The number in the
%                      flag (1 or 2) says if its P1 or P2 that is the
%                      actual end/free point. 
%
% Dn,Ds              - Displacement discontinity at each element in the
%                      model. (normal and shear). 
%
% HalfLength         - The halflength of each element in the model. 
%
%       E            - Youngs modulus.
%
%       nu           - The Poisson's ratio.
%
% LineNormalVector   - The direction cosines, CosAx (Nx) and CosAy in a list
%                      for each element.  
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
%    MidPoint      - The x and y locations of each elements midpoint. 
%
% Arguments: (output)
% K1,K2             - A vector of NaNs the length of Dn and Ds. Elements
%                     that are end elements (as in input flag) will have
%                     the value of stress intensity.
%
% Example usage:
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University

EndElsFlg=EndElsLoc>0;

%Init Array. 
K1=zeros(size(Dn));
K2=zeros(size(Dn));

ElLength=HalfLength*2;
%Eq 16 Ritz and Pollard (2012). Note its not half length that has to be
%used. 
Eq=(((0.798).*sqrt(pi).*E)./(4.*(1-(nu.^2)).*sqrt(ElLength(EndElsFlg))));

%Calculate K1 and K2. 
K1(EndElsFlg)=Dn(EndElsFlg).*Eq;
K2(EndElsFlg)=Ds(EndElsFlg).*Eq;

%% Correct sign of shear components
% This means these match the drawings in Fig 9.30 of Pollard and Fletcher
% assuming our end element normal corresponds to the y-axis in this figure.

% %Get the index of the end elements
% [Indx]=find(EndElsLoc~=0);
% for i =1:numel(Indx)
%     I=Indx(i); %Get current index
%     if EndElsLoc(Indx(i))==1
%         %Vector from midpoint to edge
%         FeM2Ev=normr([(P1(I,1)-MidPoint(I,1)),(P1(I,2))-MidPoint(I,2)]);        
%     else %must be 2
%         FeM2Ev=normr([P2(I,1)-(MidPoint(I,1)),P2(I,2)-(MidPoint(I,2))]);        
%     end
% 
%     %Rotate line normal counter clockwise 90 degrees
%     %Clockwise [x,y]=[-y,x]
%     ClockwisePerp(:,1)=-LineNormalVector(I,2);
%     ClockwisePerp(:,2)=LineNormalVector(I,1);
% 
%     %Check the two vectors allign
%     Vect=(dot(FeM2Ev',ClockwisePerp'))<=0;
%     %If they do then change sign of calculated K2 value. 
%     if Vect==1
%         K2(I,:)=-K2(I,:);
%     end
% end

%A positive shearing line dislocation in this code results in left lateral
%shear. This is the opposite of stress intensity tip displacement
%conventions
K2=-K2; 



end


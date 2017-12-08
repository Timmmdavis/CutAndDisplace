function PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )
% PlotOpeningVsShearOnEls: Function that draws the opening and shearing at 
%                         each element on the fracture surface as a quiver
%                         direction. (real coordiantes)
%
%                          If the blue vector is directed opposite to the
%                          normal direction then the element is closing.
%
%               
% usage #1: Prints results to cmd window
% PlotOpeningVsShearOnEls( LineNormalVector,Dn,Ds,P1,P2,MidPoint,HalfLength )
%
%
% Arguments: (input)
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element.  
%
%     Dn,Ds        - Two vectors representing the normal (Dn) and shear
%                   (Ds) displacement on each element.
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
%    MidPoint      - The x and y locations of each elements midpoint. 
%
%  HalfLength      - The half length of each element
%
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Grab direction cosines of each element
CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2);

%Compute the directions of the slip components.
DnCosine=LineNormalVector;
DsCosine=[-CosAy,CosAx];

%Cartesian directions and magnitudes of slip
Dn_xy=(bsxfun(@times,Dn,DnCosine));
Ds_xy=(bsxfun(@times,Ds,DsCosine));

%Summed to get total Cartesian displacement.
Dispxy=Dn_xy+Ds_xy;

% Plot slip vectors on the elements as a quiver plot. This plots each
% fracture as a line and the direction of slip for each element on the
% fracture
figure,hold on
PlotFracture( P1,P2,'r' )
xe=MidPoint(:,1);ye=MidPoint(:,2);
axis equal,
quiver(xe+(CosAx.*HalfLength),ye+(CosAy.*HalfLength),Dispxy(:,1),Dispxy(:,2),'color','b'); %normal side blue
quiver(xe-(CosAx.*HalfLength),ye-(CosAy.*HalfLength),-Dispxy(:,1),-Dispxy(:,2),'color','g'); %other side green
hold off
title('Elemental displacement directions, Element Normal Side=Blue, vectors not to scale'), xlabel('x'), ylabel('y')

%Make figure white
WhiteFigure

end


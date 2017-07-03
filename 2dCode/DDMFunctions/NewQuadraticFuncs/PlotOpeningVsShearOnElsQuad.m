function PlotOpeningVsShearOnElsQuad( NormAng,TensileDisp,ShearDisp,Points,x,y,HalfLength )
%PlotOpeningVsShearOnEls draws the directions elements slip in (world
%coords)
%As its quads Im just drawing these at the centre collation pnt

%Magnitude of shape functions
N=3;


%   Copyright 2017, Tim Davis, The University of Aberdeen
axDeg=(NormAng*180)/pi;
BDeg=axDeg+90;
u_Ds=cosd(BDeg).*ShearDisp(2:N:end,:);         %Tangent (X)
v_Ds=sind(BDeg).*ShearDisp(2:N:end,:);         %opp (y)
u_Dn=cosd(axDeg).*TensileDisp(2:N:end,:);      %Tangent (X)
v_Dn=sind(axDeg).*TensileDisp(2:N:end,:);      %opp (y)
u=u_Ds-u_Dn;
v=v_Ds-v_Dn; 
ax=NormAng;
nx=cos(ax);
ny=cos((pi/2)-ax);

% Plot slip vectors on the elements as a quiver plot. This plots each
% fracture as a line and the direction of slip for each element on the
% fracture
% Note this is only the slip one direction. 
figure,hold on
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
axis equal,
quiver(x(2:N:end,:)+(nx.*HalfLength),y(2:N:end,:)+(ny.*HalfLength),u,v,'color','b'); %normal side blue
quiver(x(2:N:end,:)-(nx.*HalfLength),y(2:N:end,:)-(ny.*HalfLength),-u,-v,'color','g'); %other side green
hold off
title('Elemental displacement directions, Element Normal Side=Blue, vectors not to scale'), xlabel('x'), ylabel('y')


end


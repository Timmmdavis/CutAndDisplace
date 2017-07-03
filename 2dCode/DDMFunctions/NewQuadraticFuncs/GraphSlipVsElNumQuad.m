function GraphSlipVsElNumQuad( NUM,TensileDisp,ShearDisp )
%GraphSlipVsElNum draws a plot of slip vs the element number. 
%The slips as shown as stair cases across the elements

x=linspace(-1,0,50);

%Magnitude of shape functions
N=3;
N2=1;
N1DnM = TensileDisp(1:N:end,:); %Normal Mag
N2DnM = TensileDisp(2:N:end,:);
N3DnM = TensileDisp(3:N:end,:);

N1DsM = ShearDisp(1:N:end,:); %Shear Mag
N2DsM = ShearDisp(2:N:end,:);
N3DsM = ShearDisp(3:N:end,:);

%Linearly increasing vec. 
Vector=zeros(size(TensileDisp));
Vector2=(1:1:(numel(TensileDisp))/3)';

Left=-(sqrt(3))/2;
Mid=0;
Right=(sqrt(3))/2;  
Left=(Left/2)-0.5;
Mid=(Mid/2)-0.5;
Right=(Right)/2-0.5;



Vector(1:N:end,:)=Vector(1:N:end,:)+Left+Vector2(1:N2:end);
Vector(2:N:end,:)=Vector(2:N:end,:)+Mid+Vector2(1:N2:end);
Vector(3:N:end,:)=Vector(3:N:end,:)+Right+Vector2(1:N2:end,:);

figure;
hold on
plot(Vector,TensileDisp,'b');
plot(Vector,ShearDisp,'r');
title('displacement discontinuity at collation pnts'), xlabel('element number')
ylabel('Dn (blue) and Ds (red) (m)'); hold off

% figure;
% plot(1:NUM*3,TensileDisp );hold on
% plot(1:NUM*3,ShearDisp );

% % Plot slip vectors on the elements as a quiver plot. This plots each
% % fracture as a line and the direction of slip for each element on the
% % fracture
% % Note this is only the slip one direction. 
% figure,hold on
% line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
% axis equal,
% quiver(x+(nx.*HalfLength),y+(ny.*HalfLength),u,v,'color','b'); %normal side blue
% quiver(x-(nx.*HalfLength),y-(ny.*HalfLength),-u,-v,'color','g'); %other side green
% hold off
% title('Elemental displacement directions, Element Normal Side=Blue, vectors not to scale'), xlabel('x'), ylabel('y')


end


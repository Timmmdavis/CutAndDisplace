function [Xnw,Ynw,Znw] = RotateObject3dNewCoords(Ax1,Ax2,Ax3,X,Y,Z)
% RotateObject3dNewCoords: Rotates the passed in XYZ points to lie in the
%               new axes defined in the inputs. Rotation assumes the centre
%               of the problem (0,0) is where we want to rotate around i.e.
%               points are already transformed correctly in space.
%               
% usage #1:
% [Xnw,Ynw,Znw] = RotateObject3dNewCoords(Ax1,Ax2,Ax3,X,Y,Z)
%
% Arguments: (input)
% X             - list of x points (mat or vect)
%
% Y             - list of y points (mat or vect)
%
% Z             - list of y points (mat or vect)
%
% Ax1           - The new orientation of the x-axis  (direction cosine 1*3)
%
% Ax2           - The new orientation of the y-axis  (direction cosine 1*3)
%
% Ax3           - The new orientation of the z-axis  (direction cosine 1*3)
%
% Arguments: (output)
% Xnw,Ynw,Znw   - the new X,Y and Z point values. 
% 
%
% Example usage:
% 
% %Ball of 1000 randomly distributed points on sphere.
% n=10000;
% r = randn(3,n); % Use a large n
% r = bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
% X = (r(1,:))';
% Y = (r(2,:))';
% Z = (r(3,:))';
% bad=Z>0;  X(bad)=[];Y(bad)=[];Z(bad)=[];
% scatter3(X,Y,Z,'b'); hold on
% %Current cosines of coord axis
% Vectx=[1,0,0];
% Vecty=[0,1,0];
% Vectz=[0,0,1];
% %Draw a RGB = XYZ axis (At 1 1 1)
% text(1,1,1,'Original coord axes') 
% quiver3(1,1,1,Vectx(1),Vectx(2),Vectx(3),'r') 
% quiver3(1,1,1,Vecty(1),Vecty(2),Vecty(3),'g')
% quiver3(1,1,1,Vectz(1),Vectz(2),Vectz(3),'b')
% %New cosine
% % Vecta=Vectx;%[1,0,0];
% % Vectb=Vecty;%[0,0,1];
% % Vectc=-Vectz;%[0,1,0];
% Vecta=[1,0,0];
% Vectb=[0,0,1];
% Vectc=[0,1,0];
% %Draw a RGB = XYZ axis NEW (At 1 2 1)
% text(1,2,1,'New coord axes') 
% quiver3(1,2,1,Vecta(1),Vecta(2),Vecta(3),'r') 
% quiver3(1,2,1,Vectb(1),Vectb(2),Vectb(3),'g')
% quiver3(1,2,1,Vectc(1),Vectc(2),Vectc(3),'b')
% [nwX,nwY,nwZ] = RotateObject3dNewCoords(Vecta,Vectb,Vectc,X,Y,Z);
% scatter3(nwX,nwY,nwZ,'r'); xlabel('x');ylabel('y');zlabel('z');
% 
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


% %Matrix multiplication with indexing 
Xnw=(Ax1(1).*X)+(Ax2(1).*Y)+(Ax3(1).*Z);
Ynw=(Ax1(2).*X)+(Ax2(2).*Y)+(Ax3(2).*Z);
Znw=(Ax1(3).*X)+(Ax2(3).*Y)+(Ax3(3).*Z);

% %Old slow way of doing this: 
% %Eq 2.23, Pollard, arranging cosines in table
% Quat=[Ax1(1),Ax2(1),Ax3(1); %x
%       Ax1(2),Ax2(2),Ax3(2); %y
%       Ax1(3),Ax2(3),Ax3(3)];%z
% %Putting these in col vecs inside var to apply the transformation 
% Points=[X(:),Y(:),Z(:)];
% %Preallocating array for loop
% NwPnts=zeros(size(Points));
% %Running loop
% for j=1:numel(X)
%     NwPnts(j,:)=(Quat*Points(j,:)')'; %Eq 2.24 Pollard and Fletcher
% end
% %Reshaping back into input array size
% Xnw=reshape(NwPnts(:,1),size(X));
% Ynw=reshape(NwPnts(:,2),size(X));
% Znw=reshape(NwPnts(:,3),size(X));
    
    
end



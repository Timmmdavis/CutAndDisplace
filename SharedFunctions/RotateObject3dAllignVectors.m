function [X,Y,Z] = RotateObject3dAllignVectors(A,B,X,Y,Z,x0,y0,z0)
% RotateObject3dAllignVectors: Alligns vector A with another vector B. 
%                 Points in the matrix are also rotated (XYZ). These are
%                 rotated around the midpoint location x0,y0,z0.
%               
% usage #1:
% [X,Y,Z] = RotateObject3dAllignVectors(A,B,X,Y,Z,x0,y0,z0)
%
% Arguments: (input)
% A             - 1x3 dir cosines, [Ax,Ay,Az]
%
% B             - 1x3 dir cosines, Ax2,Ay2,Az2 to which A is to be rotated 
%                 to allign with.
%
% X,Y,Z         - values of points that will be rotated along with the 
%                 vector A.
%
% x0,y0,z0      - rotation midpoint, single values
%
% Arguments: (output)
% X,Y,Z           - New [X,Y,Z] point locations.
% 
%
% Example usage (1):
%
% %Ball of 1000 randomly distributed points on sphere.
% n=1000;
% r = randn(3,n); % Use a large n
% r = bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
% X = (r(1,:))';
% Y = (r(2,:))';
% Z = (r(3,:))';
% scatter3(X,Y,Z,'b'); hold on
% %The point we rotate around.
% x0=0;
% y0=0;
% z0=1; 
% V1=[0,0,1 ]; %Pointing up
% V2=[1,0,0]; %Pointing East
% %Draw these
% quiver3(x0,y0,z0,V1(1),V1(2),V1(3),'b')
% quiver3(x0,y0,z0,V2(1),V2(2),V2(3),'r')
% %Now allign
% [X,Y,Z] = RotateObject3dAllignVectors(V1,V2,X,Y,Z,x0,y0,z0);
% scatter3(X,Y,Z,'r'); 
%
% Example usage (2):
%
% %Alligning triangle (normal) with Z axis, in BEM code:
% [Points(:,2),Points(:,3),Points(:,4)] = RotateObject3dAllignVectors(...
% FaceNormalVector,[0,0,1],Points(:,2),Points(:,3),Points(:,4),MidPoint(1),MidPoint(2),MidPoint(3));
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%making sure these are col vecs
X=X(:);Y=Y(:);Z=Z(:);

%midpoint at 0,0,0
X=X-x0;
Y=Y-y0;
Z=Z-z0;

%makes sure A and B are col vecs
[ A,B ] = RowVecToCol( A,B );

%The formula below doesnt work for opposing vecs so:
if B==-A
    X=-X;
    Y=-Y;
    Z=-Z;
    %midpoint back to orig Loc
    X=X+x0;
    Y=Y+y0;
    Z=Z+z0;
    return %Leave func
end

%Rotation Matrix between the two vectors
%http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
% 1. rotation vector
w=cross(A,B);
w=w/norm(w);
if any(isnan(w))
    w=[0;0;0];
end
w_hat=fcn_GetSkew(w);
% 2. rotation angle
cos_tht=A'*B/norm(A)/norm(B);
tht=acos(cos_tht);
% 3. rotation matrix, using Rodrigues' formula
r=eye(size(A,1))+w_hat*sin(tht)+w_hat^2*(1-cos(tht));


%prepping array
NEWXYZ=zeros(numel(X),3);

%rotating points
for i=1:numel(X)
    %actual rotation
    XYZ=r*[X(i);Y(i);Z(i)]; 
    %Allocating this
    NEWXYZ(i,:)=XYZ';
end

%putting back into col vecs
X=NEWXYZ(:,1);
Y=NEWXYZ(:,2);
Z=NEWXYZ(:,3);

%midpoint back to orig Loc
X=X+x0;
Y=Y+y0;
Z=Z+z0;

%internal func
function x_skew=fcn_GetSkew(x)
    x_skew=[[0,    -x(3),  x(2)]
            [x(3),     0, -x(1)]
            [-x(2), x(1),  0   ]];
end

end



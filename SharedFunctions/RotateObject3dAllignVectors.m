function [X,Y,Z] = RotateObject3dAllignVectors(A,B,X,Y,Z,x0,y0,z0)
%Alligns vector A with another vector B. Points in the matrix are also
%rotated (XYZ). These are rotated around the midpoint location x0y0z0.

%A - 1x3 dir cosines, [Ax,Ay,Az]
%B - 1x3 dir cosines, Ax2,Ay2,Az2 to which A is to be rotated onto
%X,Y,Z values of points that will be rotated along with the vector A
%x0,y0,z0 rotation midpoint, single values
%Outputs XYZ - New XYZ point locations 

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Example, alligning triangle (normal) with Z axis, code:
%[Points(:,2),Points(:,3),Points(:,4)] = RotateObject3dAllignVectors(FaceNormalVector,[0,0,1],Points(:,2),Points(:,3),Points(:,4),MidPoint(1),MidPoint(2),MidPoint(3));

%making sure these are col vecs
X=X(:);Y=Y(:);Z=Z(:);

%midpoint at 0,0,0
X=X-x0;
Y=Y-y0;
Z=Z-z0;

%makes sure A and B are col vecs
[ A,B ] = RowVecToCol( A,B );

%Rotation Matrix between the two vectors
%http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
% 1. rotation vector
w=cross(A,B);
w=w/norm(w);
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
x_skew=[0 -x(3) x(2);
 x(3) 0 -x(1);
 -x(2) x(1) 0];
end
end



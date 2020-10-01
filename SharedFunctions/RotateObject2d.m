function [Xrot,Yrot] = RotateObject2d(X,Y,Theta)
% RotateObject2d: Rotates the passed in XY points around the origin by
%               angle theta (radians) anticlockwise. Works on matricies or vectors
%               
% usage #1:
% [Xrot,Yrot] = RotateObject2d(X,Y,Theta)
%
% Arguments: (input)
% X             - list of x points (mat or vect)
%
% Y             - list of y points (mat or vect)
%
% Theta         - rotation (single value). Define in radians. 
%
% Arguments: (output)
% Xrot,Yrot     - the new X and Y values 
% 
%
% Example usage:
%
% %Create a square alligned with the axes (red). 
% x = [0 0 1 1 0]; 
% y = [0 1 1 0 0]; 
% Theta=deg2rad(45);
% fill(x,y,'r'); hold on
% quiver(-1,0,0,1); text(-1,0,'Old north')
% quiver(-0.5,0,-sin(Theta),cos(Theta)); text(-0.5,0,'New "north"')
% [Xrot,Yrot] = RotateObject2d(x,y,Theta);
% %Draw this rotated 45 degrees around the origin (blue). 
% fill(Xrot,Yrot,'b','FaceAlpha',.3); 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


Ct=cos(Theta);
St=sin(Theta);

%Just doing matrix multiplication with indexing.
Xrot=(Ct*X)+(-St*Y);
Yrot=(St*X)+(Ct*Y);


% %Old slow way of doing this: 
% %Preallocating blank arrays of right size 
% %Eq 2.23, Pollard, arranging cosines of new directions in table
% Xrot=zeros(size(X));   %Prr in polar
% Yrot=Xrot;             %Ptt in polar
% Quat=[Ct,-St;
%       St,Ct];
% 
% for i=1:numel(X)
% 
%     %Chucking the tensor together.
%     Tensor=[X(i);Y(i)];
% 
%     %http://continuummechanics.org/stressxforms.html
%     %Computing rotation
%     RotatedXY=Quat*Tensor;
% 
%     Xrot(i)=RotatedXY(1);
%     Yrot(i)=RotatedXY(2);
% 
% end

end

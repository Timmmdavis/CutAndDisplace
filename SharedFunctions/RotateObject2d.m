function [Xrot,Yrot] = RotateObject2d(X,Y,Theta )
%Rotates the passed in XY points around the origin by angle theta

%X = list of x points
%Y = list of y points
%Theta = rotation (single value). Define in radians. 
%Xrot,Yrot = the new X and Y values

%   Copyright 2017, Tim Davis, The University of Aberdeen

%making sure everything is col vecs and col vectors not grids 
%[ X,Y ] = RowVecToCol( X(:),Y(:) );

%Performing calculation for cartesian stresses
%Eq 2.23, Pollard, arranging cosines of new directions in table
Ct=cos(Theta);
St=sin(Theta);
%Quat=[Ct,-St;St,Ct];%dir 2, has to be 90 to 1 so is this. 


Xrot=(Ct*X)+(-St*Y);
Yrot=(St*X)+(Ct*Y);

%Old way of doing this: 
% % %Preallocating blank arrays of right size 
% Xrot=zeros(size(X));   %Prr in polar
% Yrot=Xrot;             %Ptt in polar

% for i=1:numel(X)
% 
% 
% %Chucking the tensor together.
% Tensor=[X(i);Y(i)];
% 
% %http://continuummechanics.org/stressxforms.html
% %Computing rotation
% RotatedXY=Quat*Tensor;
% 
% Xrot(i)=RotatedXY(1);
% Yrot(i)=RotatedXY(2);
% 
% end

end

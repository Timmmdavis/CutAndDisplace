function [RotatedCosine]=RotateCosine3d(InputCosine,Angle,RotationAxis)
% RotateCosine3d: Rotate a direction cosine/s around a specified Cartesian
%                   axis by a defined angle value (degrees). 
%               
% usage #1:
% [RotatedCosine]=RotateCosine3d(InputCosine,Angle,RotationAxis)
%
% Arguments: (input)
% InputCosine       - n*3 direction cosines, (col vec).
%
% Angle             - Angle in degress that cosines will be rotated around.
%
% RotationAxis      - A string specifying the axis. 'x','y' or 'z'
%
% Arguments: (output)
% RotatedCosine     - The input cosine but rotated. 
% 
%
% Example usage:
%
% %Ball of 1000 randomly distributed points on sphere.
% n=1000;
% r = randn(3,n); % Use a large n
% r = bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
% x = (r(1,:))';
% y = (r(2,:))';
% z = (r(3,:))';
% %Drawing this (Note its locations but also happens to be direction cosines,
% %handy!)
% quiver3(x,y,z,x,y,z); title('pointing outwards')
% %Going to rotate 90 deg around the x-axis.
% Angle=90;
% RotationAxis='x';
% [NrmVecXRot]=RotateCosine3d([x,y,z],Angle,RotationAxis);
% %Getting data out of the array
% CosAx=NrmVecXRot(:,1);
% CosAy=NrmVecXRot(:,2);
% CosAz=NrmVecXRot(:,3);
% %Drawing rotated data.
% figure;
% quiver3(x,y,z,CosAx,CosAy,CosAz); title('rotated around x-axis')
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%If we just rotate around a single angle we want to do it for all vects. 
if numel(Angle)==1
    Angle=ones(size(InputCosine(:,1)))*Angle;
end

%Preallocating arrays
RotatedCosine=zeros(size(InputCosine));
n = numel(Angle);
RotationMatrix=zeros (3,3,n);

%Filling this array with 3x3 tensors, each 3rd dim is each point
%Accumulated using indexing
if RotationAxis=='x'
        %[1,    0,      0       ]
        %[0,    cos(a), -sin(a) ]
        %[0,    sin(a), cos(a)  ]
    RotationMatrix(1,1,1:1:end) = 1;
    RotationMatrix(2,2,1:1:end) =  cos(Angle(1:1:end,:));
    RotationMatrix(2,3,1:1:end) = -sin(Angle(1:1:end,:));
    RotationMatrix(3,2,1:1:end) =  sin(Angle(1:1:end,:));
    RotationMatrix(3,3,1:1:end) =  cos(Angle(1:1:end,:));

elseif RotationAxis=='y'   
        %[cos(a),   0,  sin(a)  ]
        %[0,        1,  0       ]
        %[-sin(a),  0,  cos(a)  ]
    RotationMatrix(1,1,1:1:end) = cos(Angle(1:1:end,:));
    RotationMatrix(1,3,1:1:end) = sin(Angle(1:1:end,:));
    RotationMatrix(2,2,1:1:end) = 1;
    RotationMatrix(3,1,1:1:end) = -sin(Angle(1:1:end,:));
    RotationMatrix(3,3,1:1:end) = cos(Angle(1:1:end,:));
else % RotationAxis must =='z'    
        %[cos(a),   -sin(a),  0]
        %[sin(a),    cos(a),  0]
        %[0,             0,   1]    
    RotationMatrix(1,1,1:1:end) = cos(Angle(1:1:end,:));
    RotationMatrix(1,2,1:1:end) = -sin(Angle(1:1:end,:));
    RotationMatrix(2,1,1:1:end) = sin(Angle(1:1:end,:));
    RotationMatrix(2,2,1:1:end) = cos(Angle(1:1:end,:));
    RotationMatrix(3,3,1:1:end) = 1;
end

for i=1:numel(Angle) %its 3 cols wide  
    %Grab the current bit in the loop
    R=RotationMatrix(:,:,i);
    zp = InputCosine(i,:)';
    Rotated=R*zp;
    RotatedCosine(i,:)= Rotated';
end

end
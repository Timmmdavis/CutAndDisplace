function [RotatedCosine]=RotateCosine3d(InputCosine,Angle,RotationAxis)
RotatedCosine=zeros(size(InputCosine));

%Preallocating array
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
function [ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector,CosAx,CosAy,CosAz )
%Calculates the dipslip and strikeslip direction cosines from a column normal
%vector with 3 cols representing aX aY aZ

%FaceNormalVector n*3 vector of direction cosines. 
%CosAx,CosAy,CosAz %Split up cosines,if you do not have 'FaceNormalVector'
%but just cosine vectors you just call as: 
%[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( 1,CosAx,CosAy,CosAz )
%plotting stuff at the bottom if you want to see these visually. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

if nargin < 2
%Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
CosAx=FaceNormalVector(:,1); 
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);
end

oneslist=ones(size(FaceNormalVector(:,1)));

%Equations from Pollard Fletcher 2.102. Note these need some modification
%if to be used as real Strike/Dip Values. Fixes are given below but commented out 
AzimuthFromCosines = atan2d(CosAx,CosAy);   %extracting az (-180 to 180 degrees)
DipFromCosines = asind(-CosAz);             %extracting dip


%%%%%%%%%%%%%%%%%%%%%
%Dipslip
%Creating direction cosines for a vector pointing down the DipSlip
%direction of each triangle. First we rotate the triangles normals to face
%north, these are then rotated around the Z axis to point down dip. These
%are then rotated back to the original azimuth. 

%Rotate every vector around az so its pointing north
RotateAroundZ=degtorad(AzimuthFromCosines);
RotationAxis='z';
[PointingNorthCosine]=ApplyRotationLoop(FaceNormalVector,RotateAroundZ,RotationAxis);

%Rotate so its down dip instead of normal to the triangle (around x axis)
RotateAroundX=oneslist.*90;
RotateAroundX=degtorad(RotateAroundX);
RotationAxis='x';
[DownDipPointingNorth]=ApplyRotationLoop(PointingNorthCosine,RotateAroundX,RotationAxis);

%Rotate entire thing so its pointing along az again (Z rotation)
RotateAroundZ=degtorad(-AzimuthFromCosines);
RotationAxis='z';
[DipSlipCosine]=ApplyRotationLoop(DownDipPointingNorth,RotateAroundZ,RotationAxis);


%%%%%%%%%%%%%%
%Creating direction cosines for a vector pointing along the StrikeSlip
%direction of each triangle. This is rotated so it represents Right Hand Rule i.e.. strike
%is 90deg to the right of the dip azimuth. 

PointingAlongXAxis=zeros(size(DipSlipCosine));
PointingAlongXAxis(:,1)=-oneslist; %Each row is [1,0,0]; which is direction cosines of vect pointing along -x 
RotateAroundZ=degtorad(-AzimuthFromCosines); %done above but allocating again for saftey
RotationAxis='z';
[StrikeSlipCosine]=ApplyRotationLoop(PointingAlongXAxis,RotateAroundZ,RotationAxis);

%%%%%%%%%%%%%
%making sure flat triangles follow the same conv as the nikko TDE Script
%that is: "% Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
% as an exception, if the normal vector points upward, the strike and dip 
% vectors point Northward and Westward, whereas if the normal vector points
% downward, the strike and dip vectors point Southward and Westward, "

%Mehdi Nikkhoo's check for flat tris, we use the same conv
eZ = [0 0 1]';
for i=1:numel(FaceNormalVector(:,1))
Vstrike = cross(eZ,FaceNormalVector(i,:)');
    if norm(Vstrike)==0
    %conv as Nikkhoo
    DipSlipCosine(i,:)= [-1 0 0];      %pointing west
        if StrikeSlipCosine(i,2)>=0
            StrikeSlipCosine(i,:)= [0 1 0]; %pointing north
        else
            StrikeSlipCosine(i,:)= [0 -1 0]; %pointing south   
        end
    end
end


%Turn this on to check these are direction cosines and to drawing figures
%of dip and strike slip vector directions on the surface
% one__=(DipSlipCosine(:,1).^2)+(DipSlipCosine(:,2).^2)+(DipSlipCosine(:,3).^2);
% one__=(StrikeSlipCosine(:,1).^2)+(StrikeSlipCosine(:,2).^2)+(StrikeSlipCosine(:,3).^2);
% figure;quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3),DipSlipCosine(:,1),DipSlipCosine(:,2),DipSlipCosine(:,3))
% xlabel('x'); ylabel('y'); axis('equal'); title('DipVectorDirection');
% figure;quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3),StrikeSlipCosine(:,1),StrikeSlipCosine(:,2),StrikeSlipCosine(:,3))
% xlabel('x'); ylabel('y'); axis('equal'); title('ShearVectorDirection');

end


function [RotatedCosine]=ApplyRotationLoop(InputCosine,Angle,RotationAxis)
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

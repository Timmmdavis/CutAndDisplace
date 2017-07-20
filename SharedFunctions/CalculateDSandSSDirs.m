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
%making sure flat triangles follow the same but opp conv as the nikko TDE Script
%that is ss points NORTH and ds points WEST (when the normal points
%up). Using flipped conv as we are not using geological but engin conv. 

%flag when DS is pretty much flat cosaz=10e-14 
flag= round(DipSlipCosine(:,3),14)==0;
if any(flag) %no point doing loop unless we need to
    for i=1:numel(flag)     %loop through flag 
        if flag(i)==0       %if that one is 0 skip
            continue
        end                 %else we know this tri is flat so...
        if CosAz(i) > 0 %flat tri faces up
        DipSlipCosine(i,:)= [1 0 0];        %pointing east
        StrikeSlipCosine(i,:)= [0 -1 0];    %pointing south
        else            %tri faces down
        DipSlipCosine(i,:)= [-1 0 0];       %pointing west
        StrikeSlipCosine(i,:)= [0 1 0];     %pointing north   
        end
    end
end
%flag when DS is pretty much flat cosaz=10e-14 
% flag= round(DipSlipCosine(:,3),14)==0; %%This works with friction better
% %but not perfectly, get weird SS bits . Better to keep friction going only
% %on dipping surfaces.
% if any(flag) %no point doing loop unless we need to
%     for i=1:numel(flag)     %loop through flag 
%         if flag(i)==0       %if that one is 0 skip
%             continue
%         end                 %else we know this tri is flat so...
%         if CosAz(i) > 0 %flat tri faces up
%         DipSlipCosine(i,:)= [-1 0 0];       %pointing west
%         StrikeSlipCosine(i,:)= [0 -1 0];    %pointing south
%         else            %tri faces down
%         DipSlipCosine(i,:)= [1 0 0];        %pointing east
%         StrikeSlipCosine(i,:)= [0 1 0];     %pointing north   
%         end
%     end
% end



%Turn this on to check these are direction cosines and to drawing figures
%of dip and strike slip vector directions on the surface
% one__=(DipSlipCosine(:,1).^2)+(DipSlipCosine(:,2).^2)+(DipSlipCosine(:,3).^2);
% one__=(StrikeSlipCosine(:,1).^2)+(StrikeSlipCosine(:,2).^2)+(StrikeSlipCosine(:,3).^2);
% figure;quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3),DipSlipCosine(:,1),DipSlipCosine(:,2),DipSlipCosine(:,3))
% xlabel('x'); ylabel('y'); axis('equal'); title('DipVectorDirection');
% figure;quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3),StrikeSlipCosine(:,1),StrikeSlipCosine(:,2),StrikeSlipCosine(:,3))
% xlabel('x'); ylabel('y'); axis('equal'); title('ShearVectorDirection');

end



function [RotationMatrix]=rotx(Angle)
%Creates a rotation matrix for X axis rotation
RotationMatrix = [1,         0,                    0;
                  0,         cos(Angle), -sin(Angle);
                  0,         sin(Angle),  cos(Angle)]; %X axis Rotation
end

function [RotationMatrix]=roty(Angle)
%Creates a rotation matrix for Y axis rotation
RotationMatrix = [cos(Angle),   0, sin(Angle);
                  0,            1,          0;
                  -sin(Angle),  0, cos(Angle)]; %Y axis Rotation
end

function [RotationMatrix]=rotz(Angle)
%Creates a rotation matrix for Z axis rotation
RotationMatrix = [cos(Angle),-sin(Angle), 0;
                  sin(Angle),cos(Angle),  0;
                  0,            0,        1]; %Z axis Rotation
end

function [RotatedCosine]=ApplyRotationLoop(InputCosine,Angle,RotationAxis)
RotatedCosine=zeros(size(InputCosine));
for i=1:numel(Angle) %its 3 cols wide  
    ithAngle=Angle(i);

    if RotationAxis=='x'
    [R]=rotx(ithAngle);
    elseif RotationAxis=='y'
    [R]=roty(ithAngle);
    else % RotationAxis must =='z'
    [R]=rotz(ithAngle);
    end

    zp = InputCosine(i,:)';
    Rotated=R*zp;
    RotatedCosine(i,:)= Rotated';
end
end

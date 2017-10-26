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


%STRIKESLIP
%Creating direction cosines for a vector pointing along the StrikeSlip
%direction of each triangle. This is rotated so it represents Right Hand Rule i.e.. strike
%is 90deg to the right of the dip azimuth. 
%Just the 2D components of the normal vector (looking from above)
Res2D=[CosAx,CosAy,zeros(size(CosAx))];
%Normalise these (each row)
Res2D=normr(Res2D);
%Find vector that is 90 to this in 2D (XY)
StrikeSlipCosine=[-Res2D(:,2),Res2D(:,1),zeros(size(CosAx))];

%DIPSLIP
%Creating direction cosines for a vector pointing down the DipSlip
%direction of each triangle. 
%Calculate the new CosAz vector
DSCosAz=cos(acos(CosAz)-(pi/2));    %(90deg to original)
%Calculate the length of vector NxNy when looked above in XY
L=sin(acos(DSCosAz));     
downdip=CosAz<0; %Cos az points down so we flip the length sign
L(downdip)=-L(downdip);
%Calculate the angle a that NxNy vectors face in, this doesnt change when
%vector is rotated around dip  
a=atan2(CosAy,CosAx);   
%Calculate new lengths
DSCosAx=cos(a).*L;
DSCosAy=sin(a).*L;
DipSlipCosine=[-DSCosAx,-DSCosAy,DSCosAz];


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




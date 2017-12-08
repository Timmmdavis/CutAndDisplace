function [StrainTTotal,StrainTChg,StrainTReg]...
    =ObsPointsStressCalculator3d(Dn,Dss,Dds,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu)
% ObsPointsStressCalculator3d: Calculates the strain on points due a fault
%                              surface with known slip. This is slower than
%                              function
%                              'CalculateStressOnSurroundingPoints' but has
%                              the advantage you can store these influence
%                              matricies if you are running a loop where
%                              both the surface and obs points geometry
%                              stay the same.
%
% usage #1:
% [TotalStrainStress,StrainStressChange,StrainStressRemote,Sxx,Syy,Szz,Sxy,Sxz,Syz]...
%     =ObsPointsStressCalculator3d(Dn,Dss,Dds,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu)
%
% Arguments: (input)
%    X,Y & Z  - The observation points
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%  Dss,Dds,Dn - Vectors that describe how much the elements displace in the
%               normal (Dn) and strike slip (Dss) and dipslip (Dds)
%               directions on the elements.
%
% Sxx,Syy,Szz
% Sxy,Sxz,Syz- Remote stress if it existed on import (useful if the
%              user defined strain boundary conditions. 
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
%       nu    - The Poisson's ratio.
%
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
% Arguments: (output)
%
% StrainTTotal    - StrainTensorTotal From the strain calculation 6*n 
%                   vector, [Exx,Eyy,Ezz,Exy,Exz,Eyz]. Remote strain and
%                   strain change induced by movement of dislocations.
%
%
% StrainTChg      - StrainTensorChange From the strain calculation 6*n 
%                   vector, [Exx,Eyy,Ezz,Exy,Exz,Eyz]. Strain change
%                   induced by movement of dislocations.
%
%
% StrainTReg      - StressTensorRegionalStrain From the strain calculation 
%                   6*n vector, [Exx,Eyy,Ezz,Exy,Exz,Eyz]. Remote strain
%                   (far field tectonic)
%
% Example usage:
%
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

n = numel(X);
%Placing the defined remote driving stress/strain on arrays that match the size of the XY observation grid points 
SxxReg=zeros(n,1)+Sxx(1,1);
SyyReg=zeros(n,1)+Syy(1,1); 
SzzReg=zeros(n,1)+Szz(1,1); 
SxyReg=zeros(n,1)+Sxy(1,1); 
SxzReg=zeros(n,1)+Sxz(1,1); 
SyzReg=zeros(n,1)+Syz(1,1); 

%HookesLaw to Calcuate strain from the remote stresses that are imported into this function
[ExxReg,EyyReg,EzzReg,ExyReg,ExzReg,EyzReg] = HookesLaw3dStress2Strain( SxxReg,SyyReg,SzzReg,SxyReg,SxzReg,SyzReg,lambda,mu ) ;

%Appending into a big Cvec matrix 
StrainTReg=[ExxReg,EyyReg,EzzReg,ExyReg,ExzReg,EyzReg];

%Calculate influence matricies
[StrainInf]=CalculateInfluenceMatrices3dObsPoints(mu,lambda,X,Y,Z,P1,P2,P3,halfspace,nu);

%Accumulate arrays 
Aexx=[StrainInf.DnExx,StrainInf.DssExx,StrainInf.DdsExx]; StrainInf = rmfield(StrainInf,{'DnExx','DssExx','DdsExx'});
Aeyy=[StrainInf.DnEyy,StrainInf.DssEyy,StrainInf.DdsEyy]; StrainInf = rmfield(StrainInf,{'DnEyy','DssEyy','DdsEyy'});
Aezz=[StrainInf.DnEzz,StrainInf.DssEzz,StrainInf.DdsEzz]; StrainInf = rmfield(StrainInf,{'DnEzz','DssEzz','DdsEzz'});
Aexy=[StrainInf.DnExy,StrainInf.DssExy,StrainInf.DdsExy]; StrainInf = rmfield(StrainInf,{'DnExy','DssExy','DdsExy'});
Aexz=[StrainInf.DnExz,StrainInf.DssExz,StrainInf.DdsExz]; StrainInf = rmfield(StrainInf,{'DnExz','DssExz','DdsExz'});
Aeyz=[StrainInf.DnEyz,StrainInf.DssEyz,StrainInf.DdsEyz]; 
clear StrainInf


%Concatenate ready for equation 
A= [Aexx;Aeyy;Aezz;Aexy;Aexz;Aeyz];  clear Aexx Aeyy Aezz Aexy Aexz Aeyz

%Concatenate vector ready for equation 
B= [Dn;Dss;Dds];  

%Run linear equation system to find result of strains on the obs points
StrainTChg=A*B;
StrainTChg=reshape(StrainTChg,[],6);

%Sum to get total stress
StrainTTotal=StrainTChg+StrainTReg;

end

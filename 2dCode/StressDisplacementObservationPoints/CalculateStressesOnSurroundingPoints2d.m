function [StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
    =CalculateStressesOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,Sxx,Syy,Sxy,nu,E,halfspace...
    ,Ds,Dn,LineNormalVector)
% CalculateStressesOnSurroundingPoints2d: calculates the stress by superposition
%               This runs the calculated/defined slip for each element and
%               works out the induced stress in the surrounding points. The
%               stresses for each element are summed to a the total stress
%               for each point within the loop.
%
%               Commented lines at the base of the script allow writing of
%               a table of the stress tensors. 
%
% usage #1:
% [StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
%  =StressDispOnSurroundingPoints(x,y,MidPoint,HalfLength,Pxx,Pyy,Pxy,nu,E,halfspace...
%   ,Ds,Dn,LineNormalVector)
%
% Arguments: (input)
%      x & y  - The observation points.
%
%    MidPoint - The element midpoints in X and Y.  
%
% HalfLength  - An array of each each elements half length.
%
% Pxx,Pyy,Pxy - Remote stresses defined by user, needed to calculate total
%              stress.
%
%       nu    - The Poisson's ratio.
%
%       E     - The Young's modulus.
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%     Dn,Ds   - The calculated or defined displacement of each element. 
%
% Arguments: (output)
% StressTTotal    - StressTensorTotal. From the stress calculation 3*n 
%                   vector,[Sxx,Syy,Sxy]. Remote stress and
%                   stress change induced by movement of dislocations.
%
% StrainTTotal    - StrainTensorTotal From the strain calculation 3*n 
%                   vector, [Exx,Eyy,Exy]. Remote strain and
%                   strain change induced by movement of dislocations.
%
% StressTChg      - StressTensorChange. From the stress calculation 3*n 
%                   vector, [Sxx,Syy,Sxy]. Stress change
%                   induced by movement of dislocations.
%
% StrainTChg      - StrainTensorChange From the strain calculation 3*n 
%                   vector, [Exx,Eyy,Exy]. Strain change
%                   induced by movement of dislocations.
%
% StressTReg      - StressTensorRegionalStress. From the stress calculation
%                   3*n vector, [Sxx,Syy,Sxy]. Remote stress
%                   (far field tectonic, not gravity, this would need to be
%                   added independently).
%
% StrainTReg      - StressTensorRegionalStrain From the strain calculation 
%                   3*n vector, [Exx,Eyy,Exy]. Remote strain
%                   (far field tectonic)
%
% Example usage:
%
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

% Getting size of observation point array
n = numel(X);

%Direction cosines
CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2); 
Beta=atan2(CosAx,-CosAy);

%Placing the defined remote driving stress/strain on arrays that match the size of the XY observation grid points 
SxxReg=zeros(n,1)+Sxx(1,1);
SyyReg=zeros(n,1)+Syy(1,1); 
SxyReg=zeros(n,1)+Sxy(1,1); 
%Remote 'Strain'
[ ExxReg,EyyReg,ExyReg ] = HookesLaw2dStress2Strain( SxxReg,SyyReg,SxyReg,E,nu );

%Appending into a big Cvec matrix 
StrainTReg= [ExxReg,EyyReg,ExyReg];
StressTReg= [SxxReg,SyyReg,SxyReg];


%Setting up loop where the stresses from each element are appended to a master array
StressTChg = zeros(n,3);
 
%Setting up progress bar
if halfspace==1
    progressbar('Calculating Change in Stress on Observation XY HS') % Create figure and set starting time
else 
    progressbar('Calculating Change in Stress on Observation XY FS') % Create figure and set starting time
end
 
NumPnts=size(MidPoint(:,1),1);
%Computing effect on elements. 
if halfspace==1
    for i = 1:NumPnts
        StressTChg =StressTChg + LDstressHS(X,Y,MidPoint(i,1),MidPoint(i,2),HalfLength(i),Beta(i),Ds(i),Dn(i),nu,E);    
        progressbar(i/NumPnts) % Update figure
    end
else
    for i = 1:NumPnts
        StressTChg =StressTChg + LDstressFS(X,Y,MidPoint(i,1),MidPoint(i,2),HalfLength(i),Beta(i),Ds(i),Dn(i),nu,E); 
        progressbar(i/NumPnts) % Update figure
    end
end

%Remote 'Strain'
[ ExxChg,EyyChg,ExyChg ] = HookesLaw2dStress2Strain( StressTChg(:,1),StressTChg(:,2),StressTChg(:,3),E,nu );

StrainTChg=[ ExxChg,EyyChg,ExyChg ];

StressTTotal=StressTChg+StressTReg;
StrainTTotal=StrainTChg+StrainTReg;




end

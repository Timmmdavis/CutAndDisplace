function [StressTTotal,StressTChg,StressTReg]...
    =ObsPointsStressCalculator2d(X,Y,MidPoint,a,Pxx,Pyy,Pxy,nu,E,halfspace...
    ,Ds,Dn,LineNormalVector)
% CalculateStressOnSurroundingPoints: Calculates the stress & disp at observation
%               points by creating a influence matrix of the elements to
%               each obs point and using matrix multiplication to find the
%               stresses and displacements. Not as fast as
%               'StressDispOnSurroundingPoints.m' but if looping and the
%               problem geometry is not changing this can be a good option.
%
% usage #1:
% [StressTTotal,StressTChg,StressTReg]...
%     =ObsPointsStressCalculator2d(X,Y,MidPoint,a,Pxx,Pyy,Pxy,nu,E,halfspace...
%     ,Ds,Dn,LineNormalVector)
%
% Arguments: (input)
%      X & Y  - The observation points locations. 
%
%    MidPoint - The element midpoints in X and Y.  
%
%       a     - An array of each each elements half length
%
% Pxx,Pyy,Pxy - Remote stresses defined by user, needed to calculate total
%              stress
%
%       nu    - The Poisson's ratio
%
%       E     - The Young's modulus
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%     Dn,Ds   - The calculated or defined displacement of each element. 
%
% LineNormalVector - The direction cosines of each element [CosAx,CosAy].
%
% Arguments: (output)
% StressTTotal    - StressTensorTotal. From the stress calculation 3*n 
%                   vector,[Sxx,Syy,Sxy]. Remote stress and
%                   stress change induced by movement of dislocations.
%
%
% StressTChg      - StressTensorChange. From the stress calculation 3*n 
%                   vector, [Sxx,Syy,Sxy]. Stress change
%                   induced by movement of dislocations.
%
%
% StressTReg      - StressTensorRegionalStress. From the stress calculation
%                   3*n vector, [Sxx,Syy,Sxy]. Remote stress
%                   (far field tectonic, not gravity, this would need to be
%                   added independently).
%
%
% Example usage:
%
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



%Getting size of observation point array
 n = numel(X);

 %Creating the influence matricies
[InfMats]= CalculateInfluenceMatrices2dObsPoints(halfspace,X,Y,MidPoint,a,nu,E,LineNormalVector );
 
%Accumulate arrays 
Asxx=[InfMats.DnSxx,InfMats.DsSxx];
Asyy=[InfMats.DnSyy,InfMats.DsSyy];
Aszz=[InfMats.DnSxy,InfMats.DsSxy];
clear InfMats

%Concatenate ready for equation 
A= [Asxx;Asyy;Aszz];  
clear Asxx Asyy Aszz
 
%Concatenate vector ready for equation 
B= [Dn;Ds];  

%Run linear equation system to find result of strains on the obs points
Stress=A*B;

%Now extract variables
[ Sxx,Syy,Sxy ] = ExtractArraysFromVector( Stress );

%Extracting two arrays from the master array of induced stress
StressTChg=[Sxx,Syy,Sxy]; 
 
%Using the regional/driving stress and adding this on top of the stresses
%induced by the elements. 
StressTReg = [Pxx(1,:),Pyy(1,:),Pxy(1,:)];

StressTTotal=StressTChg+StressTReg;


end

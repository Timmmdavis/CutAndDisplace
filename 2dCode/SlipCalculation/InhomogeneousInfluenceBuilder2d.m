function [StressInfE1FB,StressInfE1IF,DispInfE1FB,DispInfE1IF,StressInfE2FB,StressInfE2IF,DispInfE2FB,DispInfE2IF,...
NUME1,NUME2,FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
= InhomogeneousInfluenceBuilder2d...
(MidPoint,HalfLength,nu,E,halfspace,LineNormalVector,BoundaryFlag)
% InhomogeneousInfluenceBuilder2d: Calculates the influence matricies for the
%               BEM calulation. This does this for both elastic parts and
%               bulds the correct flags so the later calculations know
%               which parts are free boundries, interfaces etc. 
%
% usage #1:
% See initial function call above
%
% Arguments: (input)
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%    MidPoint - The element midpoints in X and Y.
%
% HalfLength  - An array of each each elements half length.
%
%       nu    - The Poisson's ratio.
%
%       E     - The Young's modulus.
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
%   BoundaryFlag   - Flag for which part of the elastic we are dealing with,
%                   different numbers in this flag represent different bits: 
%                   0=free boundary E1
%                   1=fixed bits of the free boundary of E1
%                   2=E1-E2 interface, E1 elastic properties, normals point towards E2
%                   3=E2-E1 interface, E2 elastic properties, normals point towards E1
%                   4=If existed would be free boundary E2 
%                   5=fixed bits of the free boundary of E2
%
%
% Arguments: (output)
% %You have to understand the notation here:
%
% FB,FreeBoundaries   - Free boundary elements, elements with no additional
%                      boundary conditions
%
% IB                  - Interface boundary elements, elements with additional
%                      boundary conditions. Interface elements of Elastic 1 and 2 should have the
%                      same traction and displacement. 
%
% StressInf          -Structure containing:
% 					DsTn,DsTs,DnTn,DnTs
%					Square influence matricies of much a displacement
%                   of one element (first part of name) effects the traction
%                   on another element (last part of name).
%
% DispInf            -Structure containing:
% 					DsUx,DsUy,DnUx,DnUy
%					Square influence matricies of much a displacement
%                   of one element (first part of name) effects the
%                   displacement at the midpoint of another
%                   element (not the element displacement itself). 
%
%    Fdisp            - Elements with fixed displacements. 
%
%    NUM             - Size of influence matricies edge for one elastic
%                     part.
%
% Example usage:
%
% [ DsTn,DsTs,DnTn,DnTs,DsUx,DsUy,DnUx,DnUy]...
% = CalculateInfluenceMatrices2d(halfspace,MidPoint,HalfLength,nu,E,LineNormalVector,Fdisp )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%%%%%%%
%E1 Elastic 1
%%%%%%%
Part=1; %Which elastic we want to extract    
[MidPointE1,HalfLengthE1,nuE1,EE1,LineNormalVectorE1,NUME1,FdispE1,FB_E1,IF_E1]...
 = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector );

%E1 Elastic 1 creating big influence matrices
[ StressInfE1,DispInfE1 ] = CalculateInfluenceMatrices2d...
(halfspace,MidPointE1,HalfLengthE1,nuE1,EE1,LineNormalVectorE1,1 );
clear MidPointE1 HalfLengthE1 LineNormalVectorE1 nuE1 EE1   

%Getting the traction inf matrices for elements on the free boundary and then interface of E1
[StressInfE1FB] = ExtractData(FB_E1,1,StressInfE1);
[StressInfE1IF] = ExtractData(IF_E1,1,StressInfE1);
clear StressInfE1

%Getting the displacement inf matrices for elements on the free boundary and then interface of E1    
[DispInfE1FB] = ExtractData(FB_E1,1,DispInfE1);    
[DispInfE1IF] = ExtractData(IF_E1,1,DispInfE1);     
clear DispInfE1

%%%%%%%
%E2 Elastic 2
%%%%%%%
Part=2; %Which elastic we want to extract
[MidPointE2,HalfLengthE2,nuE2,EE2,LineNormalVectorE2,NUME2,FdispE2,FB_E2,IF_E2]...
 = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector );

%E2 Elastic 1 creating big influence matrices
[StressInfE2,DispInfE2] = CalculateInfluenceMatrices2d...
(halfspace,MidPointE2,HalfLengthE2,nuE2,EE2,LineNormalVectorE2,1 );
clear MidPointE2 HalfLengthE2 LineNormalVectorE2 nuE2 EE2 

%Getting the traction inf matrices for elements on the free boundary and then interface of E2
[StressInfE2FB] = ExtractData(FB_E2,1,StressInfE2);
[StressInfE2IF] = ExtractData(IF_E2,1,StressInfE2 );
clear StressInfE2

%Getting the displacement inf matrices for elements on the free boundary and then interface of E2    
[DispInfE2FB] = ExtractData(FB_E2,1,DispInfE2);    
[DispInfE2IF] = ExtractData(IF_E2,1,DispInfE2);     
clear DispInfE2     


%Creating some flags for fixing disps
FreeBoundaries=(BoundaryFlag == 0 | BoundaryFlag == 1 | BoundaryFlag == 4 | BoundaryFlag == 5);%Any free boundaries in either elastic   
FixedDisps=(BoundaryFlag == 1 | BoundaryFlag == 5);%Any free boundaries in either elastic 
FixedDisps=FixedDisps(FreeBoundaries); %Fixed displacements on the freeboundary elements
FreeBoundary2=(BoundaryFlag == 4 | BoundaryFlag == 5);%2nd free boundary
FreeBoundary2=FreeBoundary2(FreeBoundaries);
FreeBoundary1=~FreeBoundary2;


end



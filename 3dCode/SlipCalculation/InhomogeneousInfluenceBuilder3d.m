function [StressInfE1FB,StressInfE1IF,DispInfE1FB,DispInfE1IF,StressInfE2FB,StressInfE2IF,DispInfE2FB,DispInfE2IF,...
NUME1,NUME2,FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
= InhomogeneousInfluenceBuilder3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,BoundaryFlag)
% InhomogeneousInfluenceBuilder3d: Calculates the influence matricies for the
%               BEM calulation. This does this for both elastic parts and
%               bulds the correct flags so the later calculations know
%               which parts are free boundries, interfaces etc. 
%
% usage #1:
% See initial function call above
%
% Arguments: (input)
%
%   BoundaryFlag   - Flag for which part of the elastic we are dealing with,
%              different numbers in this flag represent different bits: 
%                   0=free boundary E1
%                   1=fixed bits of the free boundary of E1
%                   2=E1-E2 interface, E1 elastic properties, normals point towards E2
%                   3=E2-E1 interface, E2 elastic properties, normals point towards E1
%                   4=If existed would be free boundary E2 
%                   5=fixed bits of the free boundary of E2
%
%    MidPoint - The element midpoints (triangle). (col vec, XYZ)
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
% FaceNormalVector - The direction cosines, CosAx (Nx), CosAy and CosAz in 
%                   a list for each element. 
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%       nu    - The Poisson's ratio.
%
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
%                      same traction andn
%
% StressInf          -Structure containing:
%                     DnTn,DnTss,DnTds
%                     DssTn,DssTss,DssTds
%                     DdsTn,DdsTss,DdsTds
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the traction
%                     on another element (last part of name).
%
% DispInf            -Structure containing:
%                     Dn_dx,Dn_dy,Dn_dz
%                     Dss_dx,Dss_dy,Dss_dz
%                     Dds_dx,Dds_dy,Dds_dz 
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the
%                     displacement at the midpoint of another
%                     element (not the elements displacement itself). 
%
%
%       Fdisp       - Elements with fixed displacements. 
%
%        NUM        - Size of influence matricies edge for one elastic
%                     part.
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
    
if max(BoundaryFlag)<2
    error('You are modelling two elastics but the 2nd elastic has is not defined within the flag')
end
    
%%%%%%%
%E1 Elastic 1
%%%%%%%
Part=1; %Which elastic we want to extract
[MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,1 );

[ StressInfE1,DispInfE1]...
= CalculateInfluenceMatrices3d(MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,halfspace,nuE1,1);
clear MidPointE1 P1E1 P2E1 P3E1  lambdaE1   

%Getting the traction inf matrices for elements on the free boundary and then interface of E1    
[StressInfE1FB]= ExtractData(FB_E1,1,StressInfE1 );
[StressInfE1IF]= ExtractData(IF_E1,1,StressInfE1 );
clear  StressInfE1
%Getting the displacement inf matrices for elements on the free boundary and then interface of E1 
[DispInfE1FB]= ExtractData(FB_E1,1,DispInfE1);
[DispInfE1IF]= ExtractData(IF_E1,1,DispInfE1);
clear  DispInfE1


%%%%%%%
%E2 Elastic 2
%%%%%%%
Part=2; %Which elastic we want to extract
[MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,nuE2,NUME2,FdispE2,FB_E2,IF_E2]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,1 );

[StressInfE2,DispInfE2]...
= CalculateInfluenceMatrices3d(MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,halfspace,nuE2,1);
clear MidPointE2 P1E2 P2E2 P3E2  lambdaE2   

%Getting the traction inf matrices for elements on the free boundary and then interface of E2    
[StressInfE2FB]= ExtractData(FB_E2,1,StressInfE2 );
[StressInfE2IF]= ExtractData(IF_E2,1,StressInfE2 );
clear  StressInfE2
%Getting the displacement inf matrices for elements on the free boundary and then interface of E2 
[DispInfE2FB]= ExtractData(FB_E2,1,DispInfE2);
[DispInfE2IF]= ExtractData(IF_E2,1,DispInfE2);
clear  DispInfE2

%Creating some flags for fixing disps
FreeBoundaries=(BoundaryFlag == 0 | BoundaryFlag == 1 | BoundaryFlag == 4 | BoundaryFlag == 5);%Any free boundaries in either elastic   
FixedDisps=(BoundaryFlag == 1 | BoundaryFlag == 5);%Any free boundaries in either elastic 
FixedDisps=FixedDisps(FreeBoundaries); %Fixed displacements on the freeboundary elements
FreeBoundary2=(BoundaryFlag == 4 | BoundaryFlag == 5);%2nd free boundary
FreeBoundary2=FreeBoundary2(FreeBoundaries);
FreeBoundary1=~FreeBoundary2;

    
end



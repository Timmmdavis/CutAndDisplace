function [ MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,NUM,Fdisp,FB,IF  ]...
    = ExtractElasticParts3d(Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp )
% ExtractElasticParts3d: For Inhomogeneous elastics extracting
%               elements in the two different elastics using a flag.
%               Function only deals with 2 elastics currently. Used before
%               calculating the influence matricies for the elements in the
%               different elastics. A intelligent loop would need to be
%               added to get this to work with multiple elastic bodies.
%
% usage #1:
% [ MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,NUM,Fdisp,FB,IF  ]...
%     = ExtractElasticParts3d(Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp )
%
% Arguments: (input)
%    Part     - Number 1 or 2. Says if we are getting the boundaries for
%              elastic 1 or 2.
%
%BoundaryFlag - Flag for which part of the elastic we are dealing with,
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
%       nu    - The Poisson's ratio.
%
%      Fdisp  - Flag telling the user if any elements are going to be fixed
%
% Arguments: (output)
% "Same as inputs" - Any arguments the same as the inputs are simply those
%                   as described in the inputs but with the correct parts
%                   extracted for each elasic in question.
%
%       NUM        - Number of elements in the elastic that is extracted.
%
%      Fdisp       - Elements that will have fixed displacements.
%
%       FB         - Free boundary of the extracted elastic, no extra
%                   boundary conditions.
%
%       IF         - Interface elements of the extracted elastic,
%                   additional boundary conditions on these elements.
%                   
%
% Example usage:
% %Extracting elements of elastic E1:
%
%     Part=1; %Which elastic we want to extract
%     [MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
%     = ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
    
if Part==1
    
    % Flag for E1   
    Flag=BoundaryFlag<=2;
    %E1 parts
    BoundaryFlag=BoundaryFlag(Flag);
    %Flag of free boundary elements of E1 (including fixed bits )
    FB=(BoundaryFlag == 0 | BoundaryFlag == 1);
    %Flag of fixed displacement elements of E1  
    Fdisp=BoundaryFlag==1;%Fdisp on E1        
    %Flag of interface elements of E1
    IF=BoundaryFlag==2;%size of interface on E1


elseif Part==2
    
    %Flag for E2   
    Flag=BoundaryFlag>=3; 
    %E2 parts
    BoundaryFlag=BoundaryFlag(Flag);
    %Flag of interface elements of E2
    IF=BoundaryFlag==3;   
    %Flag of free boundary elements of E2 (including fixed bits )
    FB=(BoundaryFlag == 4 | BoundaryFlag == 5);
    %Flag of fixed displacement elements of E2
    Fdisp=BoundaryFlag==5; 

end

%Grabbing the elastic constants for this elastic
NUM = sum(Flag);
nu=nu(Part);
mu=mu(Part);
lambda=lambda(Part);

%Grabbing the element and thier orientations for this elastic. 
P1=P1([Flag,Flag,Flag]);     
P2=P2([Flag,Flag,Flag]);  
P3=P3([Flag,Flag,Flag]);      
MidPoint=MidPoint([Flag,Flag,Flag]);
FaceNormalVector=FaceNormalVector([Flag,Flag,Flag]);

%Reshaping to normal size now we messed with these with a flag. 
P1=reshape(P1,[],3);
P2=reshape(P2,[],3);
P3=reshape(P3,[],3);    
MidPoint=reshape(MidPoint,[],3);
FaceNormalVector=reshape(FaceNormalVector,[],3);

end


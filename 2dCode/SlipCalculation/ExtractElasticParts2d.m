function [ MidPoint,HalfLength,nu,E,LineNormalVector,NUM,Fdisp,FB,IF ]...
    = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector )
% ExtractElasticParts2d: For Inhomogeneous elastics extracting
%               elements in the two different elastics using a flag.
%               Function only deals with 2 elastics currently. Used before
%               calculating the influence matricies for the elements in the
%               different elastics. A intelligent loop would need to be
%               added to get this to work with multiple elastic bodies.
%
% usage #1:
% [ MidPoint,HalfLength,nu,E,LineNormalVector,NUM,Fdisp,FB,IF ]...
%     = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector )
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
%    MidPoint - The element midpoints in X and Y.
%
% HalfLength  - An array of each each elements half length.
%
%       nu    - The Poisson's ratio [1*2] array, nu(1) is E1 and nu(2) E2.
%
%       E     - The Young's modulus, E(1) is E1 and E(2) E2.
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
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
%     [MidPointE1,aE1,nuE1,EE1,LineNormalVectorE1,NUME1,FdispE1,FB_E1,IF_E1]...
%      = ExtractElasticParts2d( Part,BoundaryFlag,MidPoint,HalfLength,nu,E,LineNormalVector );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


    
if Part==1
    
    %Flag for E1 
    Flag=BoundaryFlag<=2;  
    %E1 parts
    BoundaryFlag=BoundaryFlag(Flag);
    %Flag of free boundary elements of E1 (including fixed bits )
    FB=(BoundaryFlag == 0 | BoundaryFlag == 1);
    %Flag of fixed displacement elements of E1  
    Fdisp=BoundaryFlag==1;   
    %Flag of interface elements of E1
    IF=BoundaryFlag==2;

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
    
%Grabbing the element orientations for this elastic. 
NUM = sum(Flag);
MidPoint=[MidPoint(Flag,1),MidPoint(Flag,2)];
HalfLength=HalfLength(Flag);  
LineNormalVector=LineNormalVector(Flag,:);  

%Grabbing the elastic constants for this elastic
nu=nu(Part);
E=E(Part);

end


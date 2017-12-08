function [StressInf,DispInf]= CalculateInfluenceMatrices2d(halfspace,MidPoint,HalfLength,nu,E,LineNormalVector,Fdisp )
% CalculateInfluenceMatrices2d: Calculates the influence matricies for the
%               BEM calulation. First this calculates the amount of stress
%               on each elements midpoint from a unit slip on each element
%               (Dn,Ds). This ends up with 2 large arrays where each of the
%               columns is the effect of one element on every midpoint of
%               the fault. When reshaped each midpoint becomes a seperate
%               row in the array. These are then converted to traction
%               influence matrices using Cauchys formula, this uses the
%               normal orientation of the fault at each midpoint.
%
% usage #1:
% [ DsTn,DsTs,DnTn,DnTs,DsUx,DsUy,DnUx,DnUy]...
% = CalculateInfluenceMatrices2d(halfspace,MidPoint,HalfLength,nu,E,LineNormalVector,Fdisp )
%
% Arguments: (input)
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%    MidPoint - The element midpoint in X and Y. 
%
%  HalfLength - An array of each each elements half length
%
%       nu    - The Poisson's ratio
%
%       E     - The Young's modulus
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
% Fdisp       - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              required).
%
%
% Arguments: (output)
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
% Example usage:
%
% [ DsTn,DsTs,DnTn,DnTs,DsUx,DsUy,DnUx,DnUy]...
% = CalculateInfluenceMatrices2d(halfspace,MidPoint,HalfLength,nu,E,LineNormalVector,Fdisp )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%If any fixed displacements exist we need to create the matrices too
FD=any(abs(Fdisp))>0;

%Number of elements
NUM=numel(MidPoint(:,1));

%Doing a memory check, will the inf matrices exceed the RAM and freeze the
%comp?
% First checking if in Octave or MATLAB, Octave has a different way of checking for free
% memory (uses Java)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    MemoryCheckerOctave(NUM,3); 
elseif isOctave==0
    if FD==0
        MemoryChecker(NUM,3); %not creating disp matrices
        %disp('function MemoryChecker is off.
    else
        MemoryChecker(NUM,5); %are creating, way bigger
        %disp('function MemoryChecker is off.
    end
end

%Get the components of the normal vector
[ CosAx,CosAy ] = ExtractCols( LineNormalVector );   

% %Moving along normal vector a very small amount to get disp correct all the
% %time.
x=MidPoint(:,1); %not moved pnts
y=MidPoint(:,2); 
x=x-(CosAx*1e-12);
y=y-(CosAy*1e-12);

%Creating empty array to be filled with influence coefficients 
Stressinfmatrix=zeros(NUM^2,3);
%Stressinfmatrix=zeros(NUM^2,5,'single');
if FD==1
    Dispinfmatrix=zeros(NUM^2,2);
    %Dispinfmatrix=zeros(NUM^2,3,'single');
else
    Dispinfmatrix=[];
end

%Setting up shear disp coeff matrix
Ds = 1;
Dn = 0;
%Running loop and filling matrices with coefficients
%Simple if/else statement to create half space or non half space
%coefficients
[DsInfMatrix,DsDisplacementXY]=CreateCoeffsLoop2d(Stressinfmatrix,...
    Dispinfmatrix,NUM,x,y,MidPoint,HalfLength,LineNormalVector,Ds,Dn,nu,E,halfspace,FD);

%Setting up normal disp coeff matrix
Ds = 0;
Dn = 1;
[DnInfMatrix,DnDisplacementXY]=CreateCoeffsLoop2d(Stressinfmatrix,...
    Dispinfmatrix,NUM,x,y,MidPoint,HalfLength,LineNormalVector,Ds,Dn,nu,E,halfspace,FD);

clear halfspace x y xe ye Beta Ds Dn first i last Stressinfmatrix Dispinfmatrix

%  Each influence array is now a huge 5*n column vectors with
%  Sxx,Syy,Szz,Ux,Uy
%  These are reshaped into 3 square matrices where each column is an
%  elements influence on every other element. 
dimx = NUM;
dimy = NUM;
%External function, extracting cols
if FD==1
    [ DsSxx,DsSyy,DsSxy ] = ExtractCols( DsInfMatrix );   
    [ DsUx,DsUy ] = ExtractCols( DsDisplacementXY ); 
    [ DsSxx,DsSyy,DsSxy,DsUx,DsUy ] = ReshapeData2d(dimx,dimy,DsSxx,DsSyy,DsSxy,DsUx,DsUy); 
else
    [ DsSxx,DsSyy,DsSxy ] = ExtractCols( DsInfMatrix );     
    [ DsSxx,DsSyy,DsSxy ] = ReshapeData2d(dimx,dimy,DsSxx,DsSyy,DsSxy);     
end
clear DsInfMatrix DsDisplacementXY

if FD==1
    [ DnSxx,DnSyy,DnSxy ] = ExtractCols( DnInfMatrix );   
    [ DnUx,DnUy ] = ExtractCols( DnDisplacementXY ); 
    [ DnSxx,DnSyy,DnSxy,DnUx,DnUy ] = ReshapeData2d(dimx,dimy,DnSxx,DnSyy,DnSxy,DnUx,DnUy); 
else
    [ DnSxx,DnSyy,DnSxy ] = ExtractCols( DnInfMatrix );     
    [ DnSxx,DnSyy,DnSxy ] = ReshapeData2d(dimx,dimy,DnSxx,DnSyy,DnSxy);     
end
clear DnInfMatrix DnDisplacementXY
clear dimx dimy

%Converts the stress influence matrices on every elements centre point to
%a XY traction influence matrix.
[ DsTx,DsTy ] = TractionVectorCartesianComponents2d( DsSxx,DsSyy,DsSxy,CosAx,CosAy );
[ DnTx,DnTy ] = TractionVectorCartesianComponents2d( DnSxx,DnSyy,DnSxy,CosAx,CosAy );
                                                    clear DsSxx DsSyy DsSxy DnSxx DnSyy DnSxy
                                                    
%Converts the traction XY to normal and shear traction components.
[ DsTn,DsTs ] = CalculateNormalShearTraction2d( DsTx,DsTy,CosAx,CosAy);
[ DnTn,DnTs ] = CalculateNormalShearTraction2d( DnTx,DnTy,CosAx,CosAy);
                                                    clear DsTx DnTx DsTy DnTy

% % %eq 5.63 C&S Element self effects, not needed but closer to real solution
%Rounding to closest self effect value (This doesn't mess up the sign convention due to movement
%along normals)
%Calculating shear mod
mu = E/(2*(1+nu));
I = eye(NUM);L = logical(I);
%C&S 5.64 element self effects
a22 = diag(HalfLength);
%Calculating self effect value
LvlstpE1=mu./((pi*(1-nu))*a22(L)); 
%Finds if intial value is positive or negative. Corrects self inf value and
%retains original sign. 
DnTn(L)=-LvlstpE1;
DsTs(L)=-LvlstpE1;
DsTn(L)=0;
DnTs(L)=0;   
                                                  
if FD==0
    [DnUx,DnUy,DsUx,DsUy ] = CreateBlankVars;        
end

%Now putting stress influence matricies inside a structure
StressInf.DnTn=DnTn;  clear DnTn
StressInf.DnTs=DnTs;  clear DnTs
StressInf.DsTn=DsTn;  clear DsTn
StressInf.DsTs=DsTs;  clear DsTs
%Now for the disp influence matricies 
DispInf.DnUx=DnUx; clear DnUx
DispInf.DnUy=DnUy; clear DnUy
DispInf.DsUx=DsUx; clear DsUx
DispInf.DsUy=DsUy; clear DsUy
                                                    
end





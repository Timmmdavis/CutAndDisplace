function [ DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy ] = CalculateInfluenceMatrices2d(NUM,halfspace,x,y,xe,ye,a,Beta,nu,E,NormAng,Fdisp )
%Calculates the influence matrices
%   First this calculates the amount of stress on each midpoint of the fault line from a unit
%   slip on each element (Sh,Ts). This ends up with 2 large arrays where
%   each of the columns is the affect of one element on every midpoint of
%   the fault. When reshaped each midpoint becomes a seperate row in the
%   array. These are then converted to traction influence matrices
%	using Cauchys formula, this uses the normal orientation of the fault at
%	each midpoint. 
%   This calls C&S functions to create list of
%   cartesian stresses induced by each element onto the other midpoints.

%   Copyright 2017, Tim Davis, The University of Aberdeen
%   x & y are element midpoints
%   xe & ye are also element midpoints, the coefficients function looks at
%   how each xe affects each x etc. 
%   a is an array of each each elements half length
%   Beta is the orientation of each element relative to the X axis in the
%   C&S convention (See C&S section 2.8 p22 for Beta definition)
%   Pxx, Pyy & Pxy are the stress inputs
%   nu is the Poisson's ratio
%   E is the Young's modulus
%   halfspace defines if we work out the coefficientsin a half or whole
%   space
%   NUM is the number of elements
%   NormAng is now the outward normal measured from Xaxis counter
%   clockwise to the normal

%If any fixed displacements exist we need to create the matrices too
FD=any(abs(Fdisp))>0;


%Doing a memory check, will the inf matrices exceed the RAM and freeze the
%comp?
% First checking if in Octave or MATLAB, Octave has a different way of checking for free
% memory (uses Java)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
%MemoryCheckerOctave2d(NUM); %doesn't work
elseif isOctave==0
MemoryChecker2d(NUM);  %Will actually be a third larger with disp matrices
end

% Computing the two direction cosines of the normal of each element of the
% fault. sin(nx) could be used instead of ny but for consistency with 3d
% code ny is used. 
ax=NormAng;
nx=cos(ax);
ny=cos((pi/2)-ax);

% %Moving along normal vector a very small amount to get disp correct all the
% %time. 
x=x-(nx*1e-12);
y=y-(ny*1e-12);

%Creating empty array to be filled with influence coefficients 
InfMatrix=zeros(NUM*NUM,5);
%InfMatrix=zeros(NUM*NUM,5,'single'); disp('Using single precision inf matrices, line 59 of CalculateInfluencematrices2d');

%Setting up shear disp coeff matrix
Ds = 1;
Dn = 0;
%Running loop and filling matrices with coefficients
%Simple if/else statement to create half space or non half space
%coefficients
StringHS='1/2 CalculatingShearDispInfMatrixHS';
StringFS='1/2 CalculatingShearDispInfMatrixFS';
[DsInfMatrix]=CreateCoeffsLoop2d(InfMatrix,...
    NUM,x,y,xe,ye,a,Beta,Ds,Dn,nu,E,halfspace,StringHS,StringFS);

%Setting up normal disp coeff matrix
Ds = 0;
Dn = 1;
StringHS='2/2 CalculatingNormalDispInfMatrixHS';
StringFS='2/2 CalculatingNormalDispInfMatrixFS';
[DnInfMatrix]=CreateCoeffsLoop2d(InfMatrix,...
    NUM,x,y,xe,ye,a,Beta,Ds,Dn,nu,E,halfspace,StringHS,StringFS);

clear halfspace x y xe ye Beta Ds Dn first i last

%  Each influence array is now a huge 5*n column vectors with
%  Sxx,Syy,Szz,ux,uy
%  Only the stress is needed
%  These are reshaped into 3 square matrices where each column is an
%  elements influence on every other element. 
dimx = NUM;
dimy = NUM;
%External function, extracting cols
if FD==1
[ DsSxx,DsSyy,DsSxy,Ds_Ux,Ds_Uy ] = ExtractCols( DsInfMatrix );     
[ DsSxx,DsSyy,DsSxy,Ds_Ux,Ds_Uy ]=ReshapeMultipleArrays(dimx,dimy,DsSxx,DsSyy,DsSxy,Ds_Ux,Ds_Uy); 
else
[ DsSxx,DsSyy,DsSxy ] = ExtractCols( DsInfMatrix );     
[ DsSxx,DsSyy,DsSxy ]=ReshapeMultipleArrays(dimx,dimy,DsSxx,DsSyy,DsSxy);     
end
clear DsInfMatrix

if FD==1
[ DnSxx,DnSyy,DnSxy,Dn_Ux,Dn_Uy ] = ExtractCols( DnInfMatrix );     
[ DnSxx,DnSyy,DnSxy,Dn_Ux,Dn_Uy ]=ReshapeMultipleArrays(dimx,dimy,DnSxx,DnSyy,DnSxy,Dn_Ux,Dn_Uy); 
else
[ DnSxx,DnSyy,DnSxy ] = ExtractCols( DnInfMatrix );     
[ DnSxx,DnSyy,DnSxy ]=ReshapeMultipleArrays(dimx,dimy,DnSxx,DnSyy,DnSxy);     
end
clear DnInfMatrix
clear dimx dimy

%Converts the stress influence matrices on every elements centre point to
%a XY traction influence matrix.
[ DsTx,DsTy ] = TractionVectorCartesianComponents2d( DsSxx,DsSyy,DsSxy,nx,ny );
[ DnTx,DnTy ] = TractionVectorCartesianComponents2d( DnSxx,DnSyy,DnSxy,nx,ny );
                                                    clear DsSxx DsSyy DsSxy DnSxx DnSyy DnSxy
                                                    
%Converts the traction XY to normal and shear traction components.
[ DsTn,DsTs ] = CalculateNormalShearTraction2d( DsTx,DsTy,nx,ny);
[ DnTn,DnTs ] = CalculateNormalShearTraction2d( DnTx,DnTy,nx,ny);
                                                    clear DsTx DnTx DsTy DnTy

% % %eq 5.63 C&S Element self effects, not needed but closer to real solution
%Rounding to closest self effect value (This doesn't mess up the sign convention due to movement
%along normals)
%Calculating shear mod
mu = E/(2*(1+nu));
I = eye(NUM);L = logical(I);
%C&S 5.64 element self effects
a22 = diag(a);
%Calculating self effect value
LvlstpE1=mu./((pi*(1-nu))*a22(L)); 
%Finds if intial value is positive or negative. Corrects self inf value and
%retains original sign. 
DnTn(L)=-LvlstpE1;
DsTs(L)=-LvlstpE1;
DsTn(L)=0;
DnTs(L)=0;   
                                                  
if FD==0
[Dn_Ux,Dn_Uy,Ds_Ux,Ds_Uy ] = CreateBlankVars;        
end
                                                    
end





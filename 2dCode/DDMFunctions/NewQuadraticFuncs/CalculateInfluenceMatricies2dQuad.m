function [ DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy ] = CalculateInfluencematrices2dQuad(NUM,halfspace,x,y,xe,ye,a,Beta,nu,E,NormAng,Fdisp )
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

%Repeating normals for quad collation pnts
[ nx3,ny3 ] = RepeatQuadValues( nx,ny );

% %Moving along normal vector a very small amount to get disp correct all the
% %time. 
x=x-(nx3*1e-12);
y=y-(ny3*1e-12);

NUM2=NUM*3;

%Creating empty array to be filled with influence coefficients 
InfMatrix=zeros(NUM2^2,5); 

G = E/(2*(1+nu));   %call outside func for efficiency

%Setting up shear disp coeff matrix
Ds = 1;
Dn = 0;
%Running loop and filling matrices with coefficients
%Simple if/else statement to create half space or non half space
%coefficients
StringHS='1/2 CalculatingShearDispInfMatrixHS';
StringFS='1/2 CalculatingShearDispInfMatrixFS';
[DsInfMatrix]=CreateCoeffsLoop(InfMatrix,...
    NUM,x,y,xe,ye,a,Beta,Ds,Dn,nu,G,halfspace,StringHS,StringFS);

%Setting up normal disp coeff matrix
Ds = 0;
Dn = 1;
StringHS='2/2 CalculatingNormalDispInfMatrixHS';
StringFS='2/2 CalculatingNormalDispInfMatrixFS';
[DnInfMatrix]=CreateCoeffsLoop(InfMatrix,...
    NUM,x,y,xe,ye,a,Beta,Ds,Dn,nu,G,halfspace,StringHS,StringFS);

clear halfspace x y xe ye Beta Ds Dn first i last

%  Each influence array is now a huge 5*n column vectors with
%  Sxx,Syy,Szz,ux,uy
%  Only the stress is needed
%  These are reshaped into 3 square matrices where each column is an
%  elements influence on every other element. 
dimx = NUM2;
dimy = NUM2;

% Now extract seperate stress column vectors into square matrices, reshaped
% so each elements influence is a seperate column. Extracting like this as
% its memory efficient
Ds_sxx = DsInfMatrix (:,1);   Ds_sxx = reshape(Ds_sxx,dimx,dimy); DsInfMatrix=DsInfMatrix(:,2:5);
Ds_syy = DsInfMatrix (:,1);   Ds_syy = reshape(Ds_syy,dimx,dimy); DsInfMatrix=DsInfMatrix(:,2:4);
Ds_sxy = DsInfMatrix (:,1);   Ds_sxy = reshape(Ds_sxy,dimx,dimy); DsInfMatrix=DsInfMatrix(:,2:3);
if FD==1
Ds_Ux  = DsInfMatrix (:,1);   Ds_Ux = reshape(Ds_Ux,dimx,dimy);   DsInfMatrix=DsInfMatrix(:,2); 
Ds_Uy  = DsInfMatrix (:,1);   Ds_Uy = reshape(Ds_Uy,dimx,dimy);   
end
clear DsInfMatrix

Dn_sxx = DnInfMatrix (:,1);   Dn_sxx = reshape(Dn_sxx,dimx,dimy); DnInfMatrix=DnInfMatrix(:,2:5);
Dn_syy = DnInfMatrix (:,1);   Dn_syy = reshape(Dn_syy,dimx,dimy); DnInfMatrix=DnInfMatrix(:,2:4);
Dn_sxy = DnInfMatrix (:,1);   Dn_sxy = reshape(Dn_sxy,dimx,dimy); DnInfMatrix=DnInfMatrix(:,2:3);
if FD==1
Dn_Ux  = DnInfMatrix (:,1);   Dn_Ux  = reshape(Dn_Ux,dimx,dimy);  DnInfMatrix=DnInfMatrix(:,2); 
Dn_Uy  = DnInfMatrix (:,1);   Dn_Uy  = reshape(Dn_Uy,dimx,dimy);
end
clear DnInfMatrix
clear dimx dimy



%Converts the stress influence matrices on every elements centre point to
%a XY traction influence matrix.
[ DsTx,DsTy ] = TractionVectorCartesianComponents2d( Ds_sxx,Ds_syy,Ds_sxy,nx3,ny3 );
[ DnTx,DnTy ] = TractionVectorCartesianComponents2d( Dn_sxx,Dn_syy,Dn_sxy,nx3,ny3 );
                                                    clear Ds_sxx Ds_syy Ds_sxy Dn_sxx Dn_syy Dn_sxy

%Converts the traction XY to normal and shear traction components.
[ DsTn,DsTs ] = CalculateNormalShearTraction2d( DsTx,DsTy,nx3,ny3);
[ DnTn,DnTs ] = CalculateNormalShearTraction2d( DnTx,DnTy,nx3,ny3);
                                                    clear ShTx ShTy DnTx DnTy
          
          
if FD==0
Dn_Ux  = [];
Dn_Uy  = [];
Ds_Ux  = [];
Ds_Uy  = [];
end


end                                                 


function [infmatrix]=CreateCoeffsLoop(infmatrix,...
    NUM,x,y,xe,ye,a,Beta,Ds,Dn,nu,G,halfspace,StringHS,StringFS)
    %Loop that calls the DDM functions and fills large column
    %coeff matrices

    %Setting up progress bar
    progressbar(StringFS)
    
    %Choosing which inf matrix to initilize
    if Ds==1
    shr=1;
    else
    shr=0;
    end
    
    %Making some vars outside the lps. 
    %Sizes for matrices
    NUM2=NUM*3;
    NUMM=NUM*9;
    %locations of the elements collation pnts (assuming el is -1 to 1). 
    Left=-(sqrt(3))/2;
    Mid=0;
    Right=(sqrt(3))/2; 
    %Shape functions for the quadratic elements. 
    %http://kratos-wiki.cimne.upc.edu/index.php/One-dimensional_Shape_Functions
    %Such that these add to 1 across whole element when used. 
    %Here numbers relate to shp funcs (n1n2n3) and letters the location on the el (coallation pnt). 
    n1lft=((-1/sqrt(3))*Left+(2/3)*Left.^2);
    n1mid=((-1/sqrt(3))*Mid+(2/3)*Mid.^2);
    n1rgt=((-1/sqrt(3))*Right+(2/3)*Right.^2);
    n2lft=(1-(4/3)*Left.^2);
    n2mid=(1-(4/3)*Mid.^2);
    n2rgt=(1-(4/3)*Right.^2);
    n3lft=((1/sqrt(3))*Left+(2/3)*Left.^2);
    n3mid=((1/sqrt(3))*Mid+(2/3)*Mid.^2);
    n3rgt=((1/sqrt(3))*Right+(2/3)*Right.^2);
    
    
    %For each element
    for i=1:NUM 
      
        %Influence matrix, how much each collation point on element i
        %effects every other element. Goes through shp func 1 2 then 3.
        InfMatrixJ=zeros(NUMM,5); 
        
        for j=1:3 %for each of the three coeff points on the els. 
        
        %Looping through the shp funcs as we create the inf matrices
        if j==1           
        mat=[n1lft,n1mid,n1rgt];
        elseif j==2
        mat=[n2lft,n2mid,n2rgt];
        elseif j==3
        mat=[n3lft,n3mid,n3rgt];
        end
        
        %Setting the Inf coeffs up
        if shr==1
        Ds=mat;
        Dn=[0,0,0];
        else
        Dn=mat;
        Ds=[0,0,0];
        end
        
        first = (j-1)*NUM2+1;
        last = j*NUM2;
                
        %Running function to create inf matrices. 
        InfMatrixJ(first:last,:)= Quad_coeff_func(x,y,xe(i),ye(i),a(i),Beta(i),Ds,Dn,nu,G);
       
        end
    
    %Adding on for each element
    first = (i-1)*NUMM+1;
    last = i*NUMM;
    infmatrix(first:last,:)=InfMatrixJ;
        
    %Updating the progress bar
    progressbar(i/NUM) 
        
    end

end



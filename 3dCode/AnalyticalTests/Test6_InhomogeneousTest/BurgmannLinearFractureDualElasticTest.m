%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%   Copyright 2017, Tim Davis, The University of Aberdeen
clear;close all

%===== add file paths ==========================
pathstring = pwd;               % Get the address of the current working directory
AdRs0=mfilename('fullpath'); %Directory containing this script
parts = strsplit(AdRs0, '\');
DirPart = parts(1,1:end-1); %Removing the file name
AdRs0 = strjoin(DirPart,'\');
cd(AdRs0) %Making sure we are in the place we are running from.

cd('..') % go up one level
cd('..') % go up one level
AdRs05 = pwd;
cd('..') % go up one level 
AdRs1 = pwd; % get the address of the new directory you're in
cd(AdRs0) % back to your original position
addpath(genpath(AdRs1))%Adds all code directorys and subfolders to the path

	
	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 


cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

	%Conventions: 
		%See figure 6.14 in David Pollards book for tensors. 
		%Geological convention is that a confining pressure is a positive stress. Geological convention is used throughout this script. 
		%For a normal stress: positive if it produces compression in the material.
		%For a shear stress: negative if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is - if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear stress Sxy will cause RightLateral movement.
		
		%Strain uses the same convention as engineering stress convention where extension is a positve value. When entering strain values remember this. 
		%For a normal stress: positive (negative) if it produces tension (compression) in the material.
		%For a shear stress: positive if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is + if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear strain Exy will cause Leftlateral movement.
		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading formatted ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %  Loads the triangles that define the fault surface 
    %  Splits this into two files that define the vertex's and which vertex's
    %  make which triangle. Use Gocad ascii format and rejig like the
    %  example files. Note column 1 for the points must start at 1.
    %  Reexport Gocad ascii file to get this. 

% S = load('TwoTriangle.dat');          %Loads the dat file -	GocadAsciiexport formatted with no text
% Points = S([1:4],[1,2,3,4]);          %Points XYZ that are the fault. First numbers define which rows, second which columns.
% Triangles = S([5:6],[1,2,3]);         %which points make a triangle
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %The file must list the traingle vertex's then triangles. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    

 
% %      %MediumSampleSetup
  string='BurgmannFreeBndE1.ts'; %Loading inner annlus, free bnd elastic 1
 [ PointsFBE1,TrianglesFBE1 ] = GoCadAsciiReader( string,pathstring );		
   string='BurgmannInterfaceE1E2.ts'; %Loading interface elastic 1
 [ PointsIFE1E2,TrianglesIFE1E2 ] = GoCadAsciiReader( string,pathstring );		
   string='BurgmannInterfaceE2E1.ts'; %Loading interface elastic 2
 [ PointsIFE2E1,TrianglesIFE2E1 ] = GoCadAsciiReader( string,pathstring );
    string='BurgmannFreeBndE2.ts'; %Loading interface elastic 2
 [ PointsFBE2,TrianglesFBE2 ] = GoCadAsciiReader( string,pathstring );


%Adding interface E1E2, Points towards E2 
[Points,Triangles,BoundaryFlag] = DataAppender3d( PointsFBE1,PointsIFE1E2,TrianglesFBE1,TrianglesIFE1E2 );

Cng=BoundaryFlag==1;
BoundaryFlag(Cng)=2; %1 represents fixed elements that do not exist in this problem

%Adding interface E2E1, Points towards E1
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsIFE2E1,Triangles,TrianglesIFE2E1,BoundaryFlag );

%Adding interface E2E1, Points towards E1
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points,PointsFBE2,Triangles,TrianglesFBE2,BoundaryFlag );


%Quick Low Down on BoundaryFlag
%0=free boundary E1
%1=fixed bits of the free boundary of E1
%2=E1-E2 interface, E1 elastic properties, normals point towards E2
%3=E2-E1 interface, E2 elastic properties, normals point towards E1
%4=If existed would be free boundary E2, 
%5=fixed bits of the free boundary of E2
%Only structured for two elastics at the moment. Introducing E3 would mean there are
%potentially 6 interfaces. E1-E2 E1-E3, E2-E1 E2=E3, E3-E1 E3-E2


 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height = 0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=100;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
lambda=mu*((2*nu)/(1-2*nu));%4000;    %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.

	%Equations you can use to calculate elastic parameters that you need above. 
	
% E = mu*(2*(1+nu)) ;               %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard

	%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
%Creating Data
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

title('All Surfaces Loaded');

	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp=(BoundaryFlag == 1 | BoundaryFlag == 5);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars3d;
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%

%Elastic 2 properties
nu2=0.25;
mu2=1000; %0.5
lambda2=mu2*((2*nu2)/(1-2*nu2));%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.

%Appending these to the back of E1
nu=[nu;nu2];
mu=[mu;mu2];
lambda=[lambda;lambda2];
  
    
Sxx = 0;       			%Positive is extension
Syy = 0; 
Szz = 0;
Sxy = 1;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;
Syz = 0; 
Option='F';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,BoundaryFlag,strain,Mu,Sf,Option);

%Removing fixed disps from BoundaryFlag
BoundaryFlag(Fdisp)=[]; %removing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
if any(Fdisp)==1
[Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)


IFFlag=(BoundaryFlag == 2 | BoundaryFlag == 3 );

ShearDispNoint=StrikeSlipDisp(~IFFlag);

PointsX=MidPoint(:,1);
PointsZ=MidPoint(:,3);
PointsXNoint=PointsX(~IFFlag);
PointsZNoint=PointsZ(~IFFlag);

MiddlePartNoInt=(PointsZNoint<100 & PointsZNoint>-100);

PointsXNointMid=PointsXNoint(MiddlePartNoInt);
ShearDispNointMid=ShearDispNoint(MiddlePartNoInt);


[Rsorted, SortIndex] = sort(PointsXNointMid);
ShearDispNointMidSort = ShearDispNointMid(SortIndex);

figure;plot(Rsorted,abs(ShearDispNointMidSort))
xlabel('x'); ylabel('y');  title('SlipDistribution'),
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Calculate distances%%%%
cells=15;   %Define the number of observations points on the grid around the fault (not extended) %6
padding=3;  %How much extra bumph you want to add away from the grid. DISTANCE not cells.        %1
[maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=MinMaxDataExtents3d(Points,cells,padding);



%%%%2d grid%%%%
[X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);    
dimx = length(X(:,1));
dimy = length(X(1,:));    
Z=zeros(size(X));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E1Flag=X<0;
XE1=X(E1Flag);
YE1=Y(E1Flag);

XE2=X(~E1Flag);
YE2=Y(~E1Flag);

%Octave can't handle transparent objects
figure;line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
title('fractures and obs points'), xlabel('x'), ylabel('y');axis equal
hold on
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
scatter(XE1(:),XE1(:),'b');
scatter(XE2(:),YE2(:),'r');
elseif isOctave==0
scatter(XE1(:),YE1(:),'b','filled','MarkerFaceAlpha',1/8);
scatter(XE2(:),YE2(:),'r','filled','MarkerFaceAlpha',1/8);
end
title('Observation Points, Blue E1 Red E2 and displacement of E1 els'), xlabel('x'), ylabel('y')
hold off         
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses  on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%TotalStress,Stress and StressChange are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
	%Strain is the first 6 coloumns and stress the last 6. Stress is the regional stress 
	%StressChange is the stress on the points from the event.
	%TotalStress is the driving stress and stress change from the event added. 
	    
%E1 Elastic 1
E1Bits =(BoundaryFlag == 0 | BoundaryFlag == 2);
TensileSlipDispE1=TensileSlipDisp(E1Bits );
StrikeSlipDispE1=StrikeSlipDisp(E1Bits );
DipSlipDispE1=DipSlipDisp(E1Bits );
Part=1; %Which elastic we want to extract
[MidPointE1,P1E1,P2E1,P3E1,muE1,lambdaE1,FaceNormalVectorE1,nuE1,NUME1,FdispE1,FB_E1,IF_E1]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [TotalStrainStressE1,StrainStressChangeE1,StrainStressRemoteE1,SyyE1,SzzE1,SxyE1,SxzE1,SyzE1]...
     =CalculateStressOnSurroundingPoints(StrikeSlipDispE1,DipSlipDispE1,TensileSlipDispE1,muE1,lambdaE1,XE1,YE1,ZE1,0,0,0,0,0,0,P1E1,P2E1,P3E1,halfspace,nuE1);

 %Extract stresses
 SxxE1=StrainStressChangeE1(:,7)';SyyE1=StrainStressChangeE1(:,8)';SxyE1=StrainStressChangeE1(:,10)';
 sz=numel(Sxx);
 
 
 %E2 Elastic 2
E2Bits =(BoundaryFlag == 3 | BoundaryFlag == 4);
TensileSlipDispE2=TensileSlipDisp(E2Bits );
StrikeSlipDispE2=StrikeSlipDisp(E2Bits );
DipSlipDispE2=DipSlipDisp(E2Bits );
Part=2; %Which elastic we want to extract
[MidPointE2,P1E2,P2E2,P3E2,muE2,lambdaE2,FaceNormalVectorE2,nuE2,NUME2,FdispE2,FB_E2,IF_E2]...
= ExtractElasticParts3d( Part,BoundaryFlag,MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,nu,Fdisp );
    
%Calculating the stress on the points that make up the annulus (elastic 1).
%0 input stress as the original analytical calc is internal pressure not
%boundary stress. its the same Tn anyway for this geometry
 [TotalStrainStressE2,StrainStressChangeE2,StrainStressRemoteE2,SyyE2,SzzE2,SxyE2,SxzE2,SyzE2]...
     =CalculateStressOnSurroundingPoints(StrikeSlipDispE2,DipSlipDispE2,TensileSlipDispE2,muE2,lambdaE2,XE2,YE2,ZE2,0,0,0,0,0,0,P1E2,P2E2,P3E2,halfspace,nuE2);

 %Extract stresses
 SxxE2=StrainStressChangeE2(:,7)';SyyE2=StrainStressChangeE2(:,8)';SxyE2=StrainStressChangeE2(:,10)';
 
